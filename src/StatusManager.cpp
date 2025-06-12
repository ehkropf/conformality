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

void StatusManager::addMessage(const StatusMessage& msg)
{
    messages.push_back(msg);

    if (messages.size() > maxMessages)
    {
        messages.erase(messages.begin());
    }
}
