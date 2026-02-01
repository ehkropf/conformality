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

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include "../src/core/StatusManager.h"

class StatusManagerTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        p_statusManager = std::make_shared<StatusManager>();
    }

    std::shared_ptr<StatusManager> p_statusManager;
};

TEST_F(StatusManagerTest, BasicMessageReporting)
{
    p_statusManager->reportDump("TestComponent", "Test dump message");
    p_statusManager->reportDebug("TestComponent", "Test debug message");
    p_statusManager->reportInfo("TestComponent", "Test info message");
    p_statusManager->reportWarning("TestComponent", "Test warning message");
    p_statusManager->reportError("TestComponent", "Test error message");

    auto messages = p_statusManager->getMessages();
    EXPECT_EQ(messages.size(), 5);

    EXPECT_EQ(messages[0].level, StatusLevel::DUMP);
    EXPECT_EQ(messages[0].component, "TestComponent");
    EXPECT_EQ(messages[0].message, "Test dump message");

    EXPECT_EQ(messages[1].level, StatusLevel::DEBUG);
    EXPECT_EQ(messages[1].component, "TestComponent");
    EXPECT_EQ(messages[1].message, "Test debug message");

    EXPECT_EQ(messages[2].level, StatusLevel::INFO);
    EXPECT_EQ(messages[2].component, "TestComponent");
    EXPECT_EQ(messages[2].message, "Test info message");

    EXPECT_EQ(messages[3].level, StatusLevel::WARNING);
    EXPECT_EQ(messages[3].component, "TestComponent");
    EXPECT_EQ(messages[3].message, "Test warning message");

    EXPECT_EQ(messages[4].level, StatusLevel::ERROR);
    EXPECT_EQ(messages[4].component, "TestComponent");
    EXPECT_EQ(messages[4].message, "Test error message");
}

TEST_F(StatusManagerTest, FilteredMessages)
{
    p_statusManager->reportDump("TestComponent", "Dump 1");
    p_statusManager->reportDebug("TestComponent", "Debug 1");
    p_statusManager->reportInfo("TestComponent", "Info 1");
    p_statusManager->reportWarning("TestComponent", "Warning 1");
    p_statusManager->reportInfo("TestComponent", "Info 2");
    p_statusManager->reportError("TestComponent", "Error 1");
    p_statusManager->reportDebug("TestComponent", "Debug 2");

    auto dumps = p_statusManager->getMessages(StatusLevel::DUMP);
    EXPECT_EQ(dumps.size(), 1);
    EXPECT_EQ(dumps[0].message, "Dump 1");

    auto debugs = p_statusManager->getMessages(StatusLevel::DEBUG);
    EXPECT_EQ(debugs.size(), 2);
    EXPECT_EQ(debugs[0].message, "Debug 1");
    EXPECT_EQ(debugs[1].message, "Debug 2");

    auto warnings = p_statusManager->getMessages(StatusLevel::WARNING);
    EXPECT_EQ(warnings.size(), 1);
    EXPECT_EQ(warnings[0].message, "Warning 1");

    auto infos = p_statusManager->getMessages(StatusLevel::INFO);
    EXPECT_EQ(infos.size(), 2);
}

TEST_F(StatusManagerTest, StatusFlags)
{
    EXPECT_FALSE(p_statusManager->hasWarnings());
    EXPECT_FALSE(p_statusManager->hasErrors());

    p_statusManager->reportDump("TestComponent", "Dump message");
    EXPECT_FALSE(p_statusManager->hasWarnings());
    EXPECT_FALSE(p_statusManager->hasErrors());

    p_statusManager->reportDebug("TestComponent", "Debug message");
    EXPECT_FALSE(p_statusManager->hasWarnings());
    EXPECT_FALSE(p_statusManager->hasErrors());

    p_statusManager->reportInfo("TestComponent", "Info message");
    EXPECT_FALSE(p_statusManager->hasWarnings());
    EXPECT_FALSE(p_statusManager->hasErrors());

    p_statusManager->reportWarning("TestComponent", "Warning message");
    EXPECT_TRUE(p_statusManager->hasWarnings());
    EXPECT_FALSE(p_statusManager->hasErrors());

    p_statusManager->reportError("TestComponent", "Error message");
    EXPECT_TRUE(p_statusManager->hasWarnings());
    EXPECT_TRUE(p_statusManager->hasErrors());
}

TEST_F(StatusManagerTest, ClearMessages)
{
    p_statusManager->reportWarning("TestComponent", "Warning message");
    p_statusManager->reportError("TestComponent", "Error message");

    EXPECT_TRUE(p_statusManager->hasWarnings());
    EXPECT_TRUE(p_statusManager->hasErrors());

    p_statusManager->clearMessages();
    EXPECT_FALSE(p_statusManager->hasWarnings());
    EXPECT_FALSE(p_statusManager->hasErrors());
    EXPECT_EQ(p_statusManager->getMessages().size(), 0);
}

TEST_F(StatusManagerTest, MaxMessagesLimit)
{
    StatusManager limitedManager(3);

    limitedManager.reportInfo("Test", "Message 1");
    limitedManager.reportInfo("Test", "Message 2");
    limitedManager.reportInfo("Test", "Message 3");
    limitedManager.reportInfo("Test", "Message 4");

    auto messages = limitedManager.getMessages();
    EXPECT_EQ(messages.size(), 3);
    EXPECT_EQ(messages[0].message, "Message 2");
    EXPECT_EQ(messages[1].message, "Message 3");
    EXPECT_EQ(messages[2].message, "Message 4");
}

TEST_F(StatusManagerTest, TieredMessageFiltering)
{
    p_statusManager->reportDump("Test", "Dump message");
    p_statusManager->reportDebug("Test", "Debug message");
    p_statusManager->reportInfo("Test", "Info message");
    p_statusManager->reportWarning("Test", "Warning message");
    p_statusManager->reportError("Test", "Error message");

    // Test filtering from DUMP level (should get all messages)
    auto fromDump = p_statusManager->getMessagesAtOrAbove(StatusLevel::DUMP);
    EXPECT_EQ(fromDump.size(), 5);

    // Test filtering from DEBUG level (should get DEBUG, INFO, WARNING, ERROR)
    auto fromDebug = p_statusManager->getMessagesAtOrAbove(StatusLevel::DEBUG);
    EXPECT_EQ(fromDebug.size(), 4);
    EXPECT_EQ(fromDebug[0].level, StatusLevel::DEBUG);
    EXPECT_EQ(fromDebug[1].level, StatusLevel::INFO);
    EXPECT_EQ(fromDebug[2].level, StatusLevel::WARNING);
    EXPECT_EQ(fromDebug[3].level, StatusLevel::ERROR);

    // Test filtering from INFO level (should get INFO, WARNING, ERROR)
    auto fromInfo = p_statusManager->getMessagesAtOrAbove(StatusLevel::INFO);
    EXPECT_EQ(fromInfo.size(), 3);
    EXPECT_EQ(fromInfo[0].level, StatusLevel::INFO);
    EXPECT_EQ(fromInfo[1].level, StatusLevel::WARNING);
    EXPECT_EQ(fromInfo[2].level, StatusLevel::ERROR);

    // Test filtering from WARNING level (should get WARNING, ERROR)
    auto fromWarning = p_statusManager->getMessagesAtOrAbove(StatusLevel::WARNING);
    EXPECT_EQ(fromWarning.size(), 2);
    EXPECT_EQ(fromWarning[0].level, StatusLevel::WARNING);
    EXPECT_EQ(fromWarning[1].level, StatusLevel::ERROR);

    // Test filtering from ERROR level (should get only ERROR)
    auto fromError = p_statusManager->getMessagesAtOrAbove(StatusLevel::ERROR);
    EXPECT_EQ(fromError.size(), 1);
    EXPECT_EQ(fromError[0].level, StatusLevel::ERROR);
}

TEST_F(StatusManagerTest, EnableLoggingRequiresFilePathForFileOutput)
{
    EXPECT_THROW(p_statusManager->enableLogging(LogOutput::FILE), std::invalid_argument);
    EXPECT_THROW(p_statusManager->enableLogging(LogOutput::BOTH), std::invalid_argument);
    EXPECT_NO_THROW(p_statusManager->enableLogging(LogOutput::NONE));
    EXPECT_NO_THROW(p_statusManager->enableLogging(LogOutput::CONSOLE));
}

TEST_F(StatusManagerTest, FileLoggingWritesToFile)
{
    std::string testLogFile = "test_status_manager.log";

    // Clean up any existing test file
    if (std::filesystem::exists(testLogFile))
    {
        std::filesystem::remove(testLogFile);
    }

    StatusManager manager;
    manager.enableLogging(LogOutput::FILE, testLogFile, StatusLevel::DEBUG);

    manager.reportInfo("TestComponent", "Test file logging message", "with details");
    manager.reportWarning("TestComponent", "Test warning");

    // Force logger to flush
    manager.flush();

    // Read the log file and verify content
    ASSERT_TRUE(std::filesystem::exists(testLogFile));

    std::ifstream logFile(testLogFile);
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    EXPECT_TRUE(content.find("TestComponent") != std::string::npos);
    EXPECT_TRUE(content.find("Test file logging message") != std::string::npos);
    EXPECT_TRUE(content.find("with details") != std::string::npos);
    EXPECT_TRUE(content.find("Test warning") != std::string::npos);

    // Clean up
    std::filesystem::remove(testLogFile);
}

TEST_F(StatusManagerTest, LogLevelFiltering)
{
    std::string testLogFile = "test_level_filter.log";

    if (std::filesystem::exists(testLogFile))
    {
        std::filesystem::remove(testLogFile);
    }

    StatusManager manager;
    // Set minimum level to WARNING (should skip DEBUG and INFO)
    manager.enableLogging(LogOutput::FILE, testLogFile, StatusLevel::WARNING);

    manager.reportDebug("Test", "Debug message");
    manager.reportInfo("Test", "Info message");
    manager.reportWarning("Test", "Warning message");
    manager.reportError("Test", "Error message");

    manager.flush();

    std::ifstream logFile(testLogFile);
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    // Debug and Info should NOT appear (filtered by level)
    EXPECT_TRUE(content.find("Debug message") == std::string::npos);
    EXPECT_TRUE(content.find("Info message") == std::string::npos);

    // Warning and Error SHOULD appear
    EXPECT_TRUE(content.find("Warning message") != std::string::npos);
    EXPECT_TRUE(content.find("Error message") != std::string::npos);

    std::filesystem::remove(testLogFile);
}

TEST_F(StatusManagerTest, MemoryStorageStillWorksWithLoggingEnabled)
{
    std::string testLogFile = "test_memory_storage.log";

    if (std::filesystem::exists(testLogFile))
    {
        std::filesystem::remove(testLogFile);
    }

    StatusManager manager;
    manager.enableLogging(LogOutput::FILE, testLogFile);

    manager.reportInfo("Test", "Message 1");
    manager.reportWarning("Test", "Message 2");
    manager.reportError("Test", "Message 3");

    // Verify messages are still stored in memory
    auto messages = manager.getMessages();
    EXPECT_EQ(messages.size(), 3);
    EXPECT_EQ(messages[0].message, "Message 1");
    EXPECT_EQ(messages[1].message, "Message 2");
    EXPECT_EQ(messages[2].message, "Message 3");

    EXPECT_TRUE(manager.hasWarnings());
    EXPECT_TRUE(manager.hasErrors());

    std::filesystem::remove(testLogFile);
}
