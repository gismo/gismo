#pragma once

class CThbGripsState
{

public:
    static CThbGripsState* GetInstance();
    static void DestroyInstance();

    GetSetPropertyMacro(ActiveLevel, int);
    GetSetPropertyMacro(ActiveRow, int);
    GetSetPropertyMacro(ActiveColumn, int);

private:
    CThbGripsState()
        : _ActiveLevel(-1)  // all
        , _ActiveRow(0)     // first row active
        , _ActiveColumn(-1) // row active 
    {
    }

    ~CThbGripsState();

    int _ActiveLevel;
    int _ActiveRow;
    int _ActiveColumn;
};

