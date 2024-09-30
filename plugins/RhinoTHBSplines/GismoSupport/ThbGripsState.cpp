#include "stdafx.h"
#include "ThbGripsState.h"

static CThbGripsState* _instance = nullptr;

CThbGripsState* CThbGripsState::GetInstance()
{
    if (!_instance)
        _instance = new CThbGripsState();
    return _instance;
}

void CThbGripsState::DestroyInstance()
{
    delete _instance;
    _instance = nullptr;
}

CThbGripsState::~CThbGripsState()
{
}
