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
#include <cassert>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace
{

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
    return spdlog::level::info;
}

} // namespace

StatusManager::StatusManager()
{
    mp_logger = spdlog::default_logger();
    assert(mp_logger != nullptr && "spdlog default logger should never be null");
}

StatusManager::StatusManager(size_t maxMsgs) : maxMessages(maxMsgs)
{
    mp_logger = spdlog::default_logger();
    assert(mp_logger != nullptr && "spdlog default logger should never be null");
}

void StatusManager::enableLogging(LogOutput output, const std::string& filePath, StatusLevel minLevel)
{
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
            sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(filePath, true));
        }
        catch (const spdlog::spdlog_ex& e)
        {
            throw std::runtime_error("Failed to create log file '" + filePath + "': " + e.what());
        }
    }

    mp_logger = std::make_shared<spdlog::logger>("conformality", sinks.begin(), sinks.end());
    mp_logger->set_level(toSpdlogLevel(minLevel));
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
    return messages;
}

std::vector<StatusMessage> StatusManager::getMessages(StatusLevel level) const
{
    std::vector<StatusMessage> filteredMessages;
    std::copy_if(messages.begin(), messages.end(), std::back_inserter(filteredMessages),
                 [level](const StatusMessage& msg) { return msg.level == level; });
    return filteredMessages;
}

std::vector<StatusMessage> StatusManager::getMessagesAtOrAbove(StatusLevel minLevel) const
{
    std::vector<StatusMessage> filteredMessages;
    std::copy_if(messages.begin(), messages.end(), std::back_inserter(filteredMessages),
                 [minLevel](const StatusMessage& msg) { return msg.level >= minLevel; });
    return filteredMessages;
}

void StatusManager::clearMessages()
{
    messages.clear();
}

bool StatusManager::hasWarnings() const
{
    return std::any_of(messages.begin(), messages.end(),
                       [](const StatusMessage& msg) { return msg.level == StatusLevel::WARNING; });
}

bool StatusManager::hasErrors() const
{
    return std::any_of(messages.begin(), messages.end(),
                       [](const StatusMessage& msg) { return msg.level == StatusLevel::ERROR; });
}

void StatusManager::flush()
{
    if (mp_logger)
    {
        mp_logger->flush();
    }
}

void StatusManager::addMessage(const StatusMessage& msg)
{
    messages.push_back(msg);

    if (messages.size() > maxMessages)
    {
        messages.erase(messages.begin());
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
}
