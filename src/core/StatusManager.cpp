/*
 * Copyright © 2025, Everett Kropf (ehkropf@gmail.com)
 *
 * This file is part of Conformality.
 * Conformality is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Conformality is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License along
 * with Conformality. If not, see <https://www.gnu.org/licenses/>.
 */

#include "StatusManager.h"

#include <algorithm>
#include <stdexcept>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace
{

// Maps StatusLevel to spdlog levels:
//   DUMP    -> trace (lowest verbosity, most detailed)
//   DEBUG   -> debug
//   INFO    -> info
//   WARNING -> warn
//   ERROR   -> err
spdlog::level::level_enum toSpdlogLevel(StatusLevel level)
{
    switch (level)
    {
        case StatusLevel::DUMP:    return spdlog::level::trace;
        case StatusLevel::DEBUG:   return spdlog::level::debug;
        case StatusLevel::INFO:    return spdlog::level::info;
        case StatusLevel::WARNING: return spdlog::level::warn;
        case StatusLevel::ERROR:   return spdlog::level::err;
    }
    // This should never be reached for valid StatusLevel values.
    // If reached, indicates memory corruption or missing enum case.
    throw std::logic_error("toSpdlogLevel: Unknown StatusLevel value: " +
                           std::to_string(static_cast<int>(level)));
}

} // namespace

StatusManager::StatusManager()
{
    // Use spdlog's default logger initially; enableLogging() can override later.
    // This ensures logging calls are safe even before explicit configuration.
    mp_logger = spdlog::default_logger();
    if (mp_logger == nullptr)
    {
        throw std::runtime_error("StatusManager: spdlog default logger is null. "
                                 "Ensure spdlog is properly initialized before creating StatusManager.");
    }
}

StatusManager::StatusManager(size_t maxMsgs) : m_maxMessages(maxMsgs)
{
    mp_logger = spdlog::default_logger();
    if (mp_logger == nullptr)
    {
        throw std::runtime_error("StatusManager: spdlog default logger is null. "
                                 "Ensure spdlog is properly initialized before creating StatusManager.");
    }
}

void StatusManager::enableLogging(LogOutput output, const std::string& filePath, StatusLevel minLevel)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_logOutput = output;

    if (output == LogOutput::NONE)
    {
        return;
    }

    std::vector<spdlog::sink_ptr> sinks;

    if (output == LogOutput::CONSOLE || output == LogOutput::BOTH)
    {
        sinks.push_back(std::make_shared<spdlog::sinks::stderr_color_sink_mt>());
    }

    if (output == LogOutput::FILE || output == LogOutput::BOTH)
    {
        if (filePath.empty())
        {
            throw std::invalid_argument("File path required for FILE or BOTH log output");
        }
        try
        {
            // Second parameter: true = truncate existing file
            sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(filePath, true));
        }
        catch (const spdlog::spdlog_ex& e)
        {
            throw std::runtime_error("Failed to create log file '" + filePath + "': " + e.what());
        }
    }

    // Defensive check: ensure at least one sink was configured
    if (sinks.empty())
    {
        throw std::logic_error("enableLogging: No sinks configured for LogOutput value: " +
                               std::to_string(static_cast<int>(output)));
    }

    mp_logger = std::make_shared<spdlog::logger>("conformality", sinks.begin(), sinks.end());
    mp_logger->set_level(toSpdlogLevel(minLevel));
    // Format: [timestamp] [colored-level] [logger-name] message
    mp_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%n] %v");
}

void StatusManager::reportDump(const std::string& component, const std::string& message,
                              const std::string& details)
{
    addMessage(StatusMessage(StatusLevel::DUMP, component, message, details));
}

void StatusManager::reportDebug(const std::string& component, const std::string& message,
                               const std::string& details)
{
    addMessage(StatusMessage(StatusLevel::DEBUG, component, message, details));
}

void StatusManager::reportInfo(const std::string& component, const std::string& message,
                              const std::string& details)
{
    addMessage(StatusMessage(StatusLevel::INFO, component, message, details));
}

void StatusManager::reportWarning(const std::string& component, const std::string& message,
                                 const std::string& details)
{
    addMessage(StatusMessage(StatusLevel::WARNING, component, message, details));
}

void StatusManager::reportError(const std::string& component, const std::string& message,
                               const std::string& details)
{
    addMessage(StatusMessage(StatusLevel::ERROR, component, message, details));
}

std::vector<StatusMessage> StatusManager::getMessages() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_messages;
}

std::vector<StatusMessage> StatusManager::getMessages(StatusLevel level) const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    std::vector<StatusMessage> filteredMessages;
    std::copy_if(m_messages.begin(), m_messages.end(), std::back_inserter(filteredMessages),
                 [level](const StatusMessage& msg) { return msg.level == level; });
    return filteredMessages;
}

std::vector<StatusMessage> StatusManager::getMessagesAtOrAbove(StatusLevel minLevel) const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    std::vector<StatusMessage> filteredMessages;
    std::copy_if(m_messages.begin(), m_messages.end(), std::back_inserter(filteredMessages),
                 [minLevel](const StatusMessage& msg) { return msg.level >= minLevel; });
    return filteredMessages;
}

void StatusManager::clearMessages()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_messages.clear();
}

bool StatusManager::hasWarnings() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return std::any_of(m_messages.begin(), m_messages.end(),
                       [](const StatusMessage& msg) { return msg.level == StatusLevel::WARNING; });
}

bool StatusManager::hasErrors() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return std::any_of(m_messages.begin(), m_messages.end(),
                       [](const StatusMessage& msg) { return msg.level == StatusLevel::ERROR; });
}

void StatusManager::flush()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (mp_logger)
    {
        mp_logger->flush();
    }
}

void StatusManager::setMaxMessages(size_t maxMsgs)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_maxMessages = maxMsgs;
}

LogOutput StatusManager::getLogOutput() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_logOutput;
}

void StatusManager::setStatusCallback(StatusCallback callback)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_statusCallback = std::move(callback);
}

void StatusManager::addMessage(const StatusMessage& msg)
{
    StatusCallback callback_copy;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_messages.push_back(msg);

        if (m_messages.size() > m_maxMessages)
        {
            m_messages.erase(m_messages.begin());
        }

        // Output via spdlog if logging is enabled
        if (m_logOutput != LogOutput::NONE && mp_logger)
        {
            std::string logMsg = "[" + msg.component + "] " + msg.message;
            if (!msg.details.empty())
            {
                logMsg += " | " + msg.details;
            }
            mp_logger->log(toSpdlogLevel(msg.level), logMsg);
        }

        callback_copy = m_statusCallback;
    }

    // Invoke callback outside the lock to avoid deadlock if callback calls back into StatusManager.
    // Wrap in try-catch so a throwing callback cannot crash the caller of reportInfo/reportWarning/etc.
    if (callback_copy)
    {
        try
        {
            callback_copy(msg);
        }
        catch (const std::exception& e)
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_logOutput != LogOutput::NONE && mp_logger)
            {
                mp_logger->error("StatusManager callback threw: {}", e.what());
            }
        }
    }
}
