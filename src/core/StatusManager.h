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

#ifndef STATUS_MANAGER_H
#define STATUS_MANAGER_H

#include <memory>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>

/**
 * @brief Configuration for log output destinations
 */
enum class LogOutput
{
    NONE,    /**< No logging output (memory storage only) */
    CONSOLE, /**< Log to console (stderr) */
    FILE,    /**< Log to file */
    BOTH     /**< Log to both console and file */
};

/**
 * @brief Enumeration of status message levels
 */
enum class StatusLevel
{
    DUMP,    /**< Debug dump information */
    DEBUG,   /**< Debug information */
    INFO,    /**< General information */
    WARNING, /**< Warning messages */
    ERROR    /**< Error messages */
};

/**
 * @brief Structure representing a status message
 */
struct StatusMessage
{
    StatusLevel level;      /**< Severity level of the message */
    std::string component;  /**< Component that generated the message */
    std::string message;    /**< Main message text */
    std::string details;    /**< Optional additional details */

    /**
     * @brief Construct a new Status Message
     * @param lvl Message severity level
     * @param comp Component name that generated the message
     * @param msg Main message text
     * @param det Optional additional details
     */
    StatusMessage(StatusLevel lvl, const std::string& comp,
                  const std::string& msg, const std::string& det = "")
        : level(lvl), component(comp), message(msg), details(det) {}
};

/**
 * @brief Abstract interface for status message management
 *
 * This interface defines the contract for status managers that handle
 * reporting and retrieval of status messages from various components.
 * It follows the dependency injection pattern for clean architecture.
 */
class IStatusManager
{
public:
    virtual ~IStatusManager() = default;

    /**
     * @brief Report a dump message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    virtual void reportDump(const std::string& component, const std::string& message,
                           const std::string& details = "") = 0;

    /**
     * @brief Report a debug message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    virtual void reportDebug(const std::string& component, const std::string& message,
                            const std::string& details = "") = 0;

    /**
     * @brief Report an informational message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    virtual void reportInfo(const std::string& component, const std::string& message,
                           const std::string& details = "") = 0;

    /**
     * @brief Report a warning message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    virtual void reportWarning(const std::string& component, const std::string& message,
                              const std::string& details = "") = 0;

    /**
     * @brief Report an error message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    virtual void reportError(const std::string& component, const std::string& message,
                            const std::string& details = "") = 0;

    /**
     * @brief Get all stored messages
     * @return Vector of all status messages
     */
    virtual std::vector<StatusMessage> getMessages() const = 0;

    /**
     * @brief Get messages filtered by level
     * @param level Status level to filter by
     * @return Vector of messages matching the specified level
     */
    virtual std::vector<StatusMessage> getMessages(StatusLevel level) const = 0;

    /**
     * @brief Get messages at or above the specified minimum level
     * @param minLevel Minimum status level (inclusive)
     * @return Vector of messages at or above the specified level
     */
    virtual std::vector<StatusMessage> getMessagesAtOrAbove(StatusLevel minLevel) const = 0;

    /**
     * @brief Clear all stored messages
     */
    virtual void clearMessages() = 0;

    /**
     * @brief Check if any warning messages exist
     * @return True if warnings are present, false otherwise
     */
    virtual bool hasWarnings() const = 0;

    /**
     * @brief Check if any error messages exist
     * @return True if errors are present, false otherwise
     */
    virtual bool hasErrors() const = 0;
};

/**
 * @brief Concrete implementation of status message management
 *
 * This class provides a concrete implementation of the IStatusManager interface,
 * managing status messages with configurable limits and filtering capabilities.
 * Messages are stored in memory with automatic cleanup when limits are exceeded.
 */
class StatusManager : public IStatusManager
{
private:
    std::vector<StatusMessage> m_messages;  /**< Storage for status messages */
    size_t m_maxMessages{1000};             /**< Maximum number of messages to store */
    std::shared_ptr<spdlog::logger> mp_logger; /**< spdlog logger instance; defaults to spdlog's default logger until enableLogging() is called */
    LogOutput m_logOutput{LogOutput::NONE};    /**< Current log output configuration */

public:
    /**
     * @brief Default constructor with default message limit
     */
    StatusManager();

    /**
     * @brief Constructor with custom message limit
     * @param maxMsgs Maximum number of messages to store
     */
    explicit StatusManager(size_t maxMsgs);

    /**
     * @brief Enable logging output
     * @param output Log output destination (CONSOLE, FILE, or BOTH)
     * @param filePath Path to log file (required if output includes FILE)
     * @param minLevel Minimum level to log (default: DEBUG, skipping DUMP)
     */
    void enableLogging(LogOutput output, const std::string& filePath = "",
                       StatusLevel minLevel = StatusLevel::DEBUG);

    /**
     * @brief Report a dump message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    void reportDump(const std::string& component, const std::string& message,
                   const std::string& details = "") override;

    /**
     * @brief Report a debug message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    void reportDebug(const std::string& component, const std::string& message,
                    const std::string& details = "") override;

    /**
     * @brief Report an informational message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    void reportInfo(const std::string& component, const std::string& message,
                   const std::string& details = "") override;

    /**
     * @brief Report a warning message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    void reportWarning(const std::string& component, const std::string& message,
                      const std::string& details = "") override;

    /**
     * @brief Report an error message
     * @param component Name of the component reporting the message
     * @param message Main message text
     * @param details Optional additional details
     */
    void reportError(const std::string& component, const std::string& message,
                    const std::string& details = "") override;

    /**
     * @brief Get all stored messages
     * @return Vector of all status messages in chronological order
     */
    std::vector<StatusMessage> getMessages() const override;

    /**
     * @brief Get messages filtered by level
     * @param level Status level to filter by
     * @return Vector of messages matching the specified level
     */
    std::vector<StatusMessage> getMessages(StatusLevel level) const override;

    /**
     * @brief Get messages at or above the specified minimum level
     * @param minLevel Minimum status level (inclusive)
     * @return Vector of messages at or above the specified level
     */
    std::vector<StatusMessage> getMessagesAtOrAbove(StatusLevel minLevel) const override;

    /**
     * @brief Clear all stored messages
     */
    void clearMessages() override;

    /**
     * @brief Check if any warning messages exist
     * @return True if warnings are present, false otherwise
     */
    bool hasWarnings() const override;

    /**
     * @brief Check if any error messages exist
     * @return True if errors are present, false otherwise
     */
    bool hasErrors() const override;

    /**
     * @brief Set the maximum number of messages to store
     * @param maxMsgs New maximum message count
     */
    void setMaxMessages(size_t maxMsgs) { m_maxMessages = maxMsgs; }

    /**
     * @brief Get the current log output configuration
     * @return Current LogOutput setting
     */
    LogOutput getLogOutput() const { return m_logOutput; }

    /**
     * @brief Flush any pending log output
     *
     * Forces any buffered log messages to be written to their destinations.
     */
    void flush();

private:
    /**
     * @brief Add a message to storage with automatic cleanup
     * @param msg Status message to add
     *
     * If the message count exceeds maxMessages, the oldest message is removed.
     */
    void addMessage(const StatusMessage& msg);
};

#endif // STATUS_MANAGER_H
