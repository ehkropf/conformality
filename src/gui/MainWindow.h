#pragma once

class MainWindow
{
private:
    bool m_showDemoWindow;
    
public:
    MainWindow();
    ~MainWindow();
    
    bool initialize();
    void render();
    void shutdown();
    
private:
    void renderMenuBar();
    void renderMainLayout();
};