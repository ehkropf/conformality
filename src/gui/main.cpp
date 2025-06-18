#include "Application.h"
#include <iostream>

int main(int argc, char** argv)
{
    Application app;
    
    if (!app.initialize())
    {
        std::cerr << "Failed to initialize application" << std::endl;
        return 1;
    }
    
    app.run();
    app.shutdown();
    
    return 0;
}