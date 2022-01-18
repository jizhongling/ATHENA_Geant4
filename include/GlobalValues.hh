/*
* these are global values that are shared between the parent and child threads
*/
#include "globals.hh"
#include <vector>
#include <tuple>

#ifndef GlobalValues_h
#define GlobalValues_h 1

namespace GlobalValues
{
    extern const G4int NumHCalLayers; // Number of total layers in HCal
    extern const G4int NumHCalTowers; // One-dimensional number of towers. Current is 6x6, so this = 6
    extern const G4int NumECalBlocks; // One-dimensional number of blocks. Current is 8x8, so this = 8
}
#endif