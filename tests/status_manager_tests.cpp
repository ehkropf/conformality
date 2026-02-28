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
#include <atomic>
#include <fstream>
#include <filesystem>
#include <thread>
#include <vector>
#include "../src/core/StatusManager.h"

// RAII helper for test file cleanup
class TestFileCleanup
{
public:
    explicit TestFileCleanup(std::string path) : m_path(std::move(path))
    {
        // Clean up any existing file at construction
        if (std::filesystem::exists(m_path))
        {
            std::filesystem::remove(m_path);
        }
    }

    ~TestFileCleanup()
    {
        if (std::filesystem::exists(m_path))
        {
            std::filesystem::remove(m_path);
        }
    }

    const std::string& path() const { return m_path; }

    // Non-copyable, non-movable
    TestFileCleanup(const TestFileCleanup&) = delete;
    TestFileCleanup& operator=(const TestFileCleanup&) = delete;

private:
    std::string m_path;
};

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
    TestFileCleanup cleanup("test_status_manager.log");

    StatusManager manager;
    manager.enableLogging(LogOutput::FILE, cleanup.path(), StatusLevel::DEBUG);

    manager.reportInfo("TestComponent", "Test file logging message", "with details");
    manager.reportWarning("TestComponent", "Test warning");

    manager.flush();

    ASSERT_TRUE(std::filesystem::exists(cleanup.path()));

    std::ifstream logFile(cleanup.path());
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    EXPECT_TRUE(content.find("TestComponent") != std::string::npos);
    EXPECT_TRUE(content.find("Test file logging message") != std::string::npos);
    EXPECT_TRUE(content.find("with details") != std::string::npos);
    EXPECT_TRUE(content.find("Test warning") != std::string::npos);
}

TEST_F(StatusManagerTest, LogLevelFiltering)
{
    TestFileCleanup cleanup("test_level_filter.log");

    StatusManager manager;
    // Set minimum level to WARNING (should skip DEBUG and INFO)
    manager.enableLogging(LogOutput::FILE, cleanup.path(), StatusLevel::WARNING);

    manager.reportDebug("Test", "Debug message");
    manager.reportInfo("Test", "Info message");
    manager.reportWarning("Test", "Warning message");
    manager.reportError("Test", "Error message");

    manager.flush();

    std::ifstream logFile(cleanup.path());
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    // Debug and Info should NOT appear (filtered by level)
    EXPECT_TRUE(content.find("Debug message") == std::string::npos);
    EXPECT_TRUE(content.find("Info message") == std::string::npos);

    // Warning and Error SHOULD appear
    EXPECT_TRUE(content.find("Warning message") != std::string::npos);
    EXPECT_TRUE(content.find("Error message") != std::string::npos);
}

TEST_F(StatusManagerTest, MemoryStorageStillWorksWithLoggingEnabled)
{
    TestFileCleanup cleanup("test_memory_storage.log");

    StatusManager manager;
    manager.enableLogging(LogOutput::FILE, cleanup.path());

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
}

TEST_F(StatusManagerTest, DefaultMinLevelFiltersDumpMessages)
{
    TestFileCleanup cleanup("test_dump_filter.log");

    StatusManager manager;
    // Enable logging with default minLevel (DEBUG) - DUMP should be filtered
    manager.enableLogging(LogOutput::FILE, cleanup.path());

    manager.reportDump("Test", "Dump should not appear");
    manager.reportDebug("Test", "Debug should appear");
    manager.reportInfo("Test", "Info should appear");

    manager.flush();

    std::ifstream logFile(cleanup.path());
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    // DUMP should NOT appear (filtered by default minLevel=DEBUG)
    EXPECT_TRUE(content.find("Dump should not appear") == std::string::npos);

    // DEBUG and INFO SHOULD appear
    EXPECT_TRUE(content.find("Debug should appear") != std::string::npos);
    EXPECT_TRUE(content.find("Info should appear") != std::string::npos);
}

TEST_F(StatusManagerTest, BothOutputWritesToFile)
{
    TestFileCleanup cleanup("test_both_output.log");

    StatusManager manager;
    manager.enableLogging(LogOutput::BOTH, cleanup.path(), StatusLevel::DEBUG);

    manager.reportInfo("Test", "Both output message");
    manager.reportWarning("Test", "Both output warning");

    manager.flush();

    // Verify file output works (console output is harder to capture in tests)
    ASSERT_TRUE(std::filesystem::exists(cleanup.path()));

    std::ifstream logFile(cleanup.path());
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    EXPECT_TRUE(content.find("Both output message") != std::string::npos);
    EXPECT_TRUE(content.find("Both output warning") != std::string::npos);
}

TEST_F(StatusManagerTest, InvalidFilePathThrowsRuntimeError)
{
    StatusManager manager;
    // Attempting to write to a non-existent directory should throw
    EXPECT_THROW(
        manager.enableLogging(LogOutput::FILE, "/nonexistent/directory/test.log"),
        std::runtime_error);
}

TEST_F(StatusManagerTest, ReconfigureLoggingOutput)
{
    TestFileCleanup cleanup1("test_reconfig_1.log");
    TestFileCleanup cleanup2("test_reconfig_2.log");

    StatusManager manager;

    // Enable FILE logging to first file
    manager.enableLogging(LogOutput::FILE, cleanup1.path());
    manager.reportInfo("Test", "Message to file 1");
    manager.flush();

    // Reconfigure to different file
    manager.enableLogging(LogOutput::FILE, cleanup2.path());
    manager.reportInfo("Test", "Message to file 2");
    manager.flush();

    // Verify Message 1 in file1
    std::ifstream logFile1(cleanup1.path());
    std::string content1((std::istreambuf_iterator<char>(logFile1)),
                         std::istreambuf_iterator<char>());
    EXPECT_TRUE(content1.find("Message to file 1") != std::string::npos);
    EXPECT_TRUE(content1.find("Message to file 2") == std::string::npos);

    // Verify Message 2 in file2
    std::ifstream logFile2(cleanup2.path());
    std::string content2((std::istreambuf_iterator<char>(logFile2)),
                         std::istreambuf_iterator<char>());
    EXPECT_TRUE(content2.find("Message to file 2") != std::string::npos);
    EXPECT_TRUE(content2.find("Message to file 1") == std::string::npos);
}

TEST_F(StatusManagerTest, FlushSafeWhenLoggingDisabled)
{
    StatusManager manager;
    // No logging enabled - flush should be safe (no-op)
    EXPECT_NO_THROW(manager.flush());
}

TEST_F(StatusManagerTest, DisableLoggingAfterEnabling)
{
    TestFileCleanup cleanup("test_disable_logging.log");

    StatusManager manager;
    manager.enableLogging(LogOutput::FILE, cleanup.path());
    manager.reportInfo("Test", "Should be logged");
    manager.flush();

    // Disable logging
    manager.enableLogging(LogOutput::NONE);
    manager.reportInfo("Test", "Should NOT be logged");
    manager.flush();

    std::ifstream logFile(cleanup.path());
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    EXPECT_TRUE(content.find("Should be logged") != std::string::npos);
    EXPECT_TRUE(content.find("Should NOT be logged") == std::string::npos);

    // Verify getLogOutput reflects the change
    EXPECT_EQ(manager.getLogOutput(), LogOutput::NONE);
}

TEST_F(StatusManagerTest, DumpLevelCanBeExplicitlyEnabled)
{
    TestFileCleanup cleanup("test_dump_enabled.log");

    StatusManager manager;
    manager.enableLogging(LogOutput::FILE, cleanup.path(), StatusLevel::DUMP);

    manager.reportDump("Test", "Dump message should appear");
    manager.flush();

    std::ifstream logFile(cleanup.path());
    std::string content((std::istreambuf_iterator<char>(logFile)),
                        std::istreambuf_iterator<char>());

    EXPECT_TRUE(content.find("Dump message should appear") != std::string::npos);
}

TEST_F(StatusManagerTest, ConsoleOnlyDoesNotCreateFile)
{
    const std::string testPath = "test_console_no_file.log";
    // Ensure file doesn't exist
    if (std::filesystem::exists(testPath))
    {
        std::filesystem::remove(testPath);
    }

    StatusManager manager;
    // Pass a file path but use CONSOLE output - file should be ignored
    manager.enableLogging(LogOutput::CONSOLE, testPath);
    manager.reportInfo("Test", "Console only message");
    manager.flush();

    // File should NOT be created for CONSOLE-only output
    EXPECT_FALSE(std::filesystem::exists(testPath));

    // Verify getLogOutput returns correct value
    EXPECT_EQ(manager.getLogOutput(), LogOutput::CONSOLE);
}

// --- Thread safety tests ---

TEST_F(StatusManagerTest, ThreadSafety_ConcurrentWrites)
{
    const int num_threads = 4;
    const int messages_per_thread = 100;

    std::vector<std::thread> threads;
    for (int t = 0; t < num_threads; ++t)
    {
        threads.emplace_back([this, t, messages_per_thread]() {
            for (int i = 0; i < messages_per_thread; ++i)
            {
                p_statusManager->reportInfo("Thread" + std::to_string(t),
                                            "Message " + std::to_string(i));
            }
        });
    }

    for (auto& thread : threads) thread.join();

    auto messages = p_statusManager->getMessages();
    EXPECT_EQ(messages.size(), static_cast<size_t>(num_threads * messages_per_thread));
}

TEST_F(StatusManagerTest, ThreadSafety_ReadWhileWriting)
{
    std::atomic<bool> stop{false};

    // Writer thread
    std::thread writer([this, &stop]() {
        int i = 0;
        while (!stop.load())
        {
            p_statusManager->reportInfo("Writer", "Message " + std::to_string(i++));
        }
    });

    // Reader thread — exercise all read paths concurrently
    std::thread reader([this, &stop]() {
        while (!stop.load())
        {
            auto all = p_statusManager->getMessages();
            auto filtered = p_statusManager->getMessages(StatusLevel::INFO);
            auto tiered = p_statusManager->getMessagesAtOrAbove(StatusLevel::WARNING);
            (void)p_statusManager->hasWarnings();
            (void)p_statusManager->hasErrors();
            (void)all;
            (void)filtered;
            (void)tiered;
        }
    });

    // Let them run for a short burst
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    stop.store(true);

    writer.join();
    reader.join();

    // No crash means success; verify we got some messages
    EXPECT_GT(p_statusManager->getMessages().size(), 0u);
}

// --- Callback tests ---

TEST_F(StatusManagerTest, Callback_InvokedOnReport)
{
    std::vector<StatusMessage> received;
    p_statusManager->setStatusCallback([&received](const StatusMessage& msg) {
        received.push_back(msg);
    });

    p_statusManager->reportInfo("TestComp", "Hello");
    p_statusManager->reportWarning("TestComp", "Warning");

    ASSERT_EQ(received.size(), 2u);
    EXPECT_EQ(received[0].level, StatusLevel::INFO);
    EXPECT_EQ(received[0].message, "Hello");
    EXPECT_EQ(received[1].level, StatusLevel::WARNING);
    EXPECT_EQ(received[1].message, "Warning");
}

TEST_F(StatusManagerTest, Callback_ClearedSafely)
{
    int callCount = 0;
    p_statusManager->setStatusCallback([&callCount](const StatusMessage&) {
        ++callCount;
    });

    p_statusManager->reportInfo("Test", "Before clear");
    EXPECT_EQ(callCount, 1);

    // Clear callback
    p_statusManager->setStatusCallback(nullptr);

    // Reporting after clearing should not crash and not invoke old callback
    EXPECT_NO_THROW(p_statusManager->reportInfo("Test", "After clear"));
    EXPECT_EQ(callCount, 1);
}
