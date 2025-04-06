#include <catch2/catch_test_macros.hpp>
#include "../include/fast_transfer_entropy.hpp"
#include <cstdlib>
#include <cmath>

TEST_CASE("Test Transfer Entropy of Identical Sequences"){
    using namespace fast_entropy;
    int numSymbols = 4;
    int maxExpectedNodes = 1000;

    SymbolSequence symbol = SymbolSequence(100);
    for (int i=0; i < 100; i++){
        symbol[i] = 0;
    }

    for (int delta = 2; delta < 10; ++delta)
    REQUIRE(calculateSymbolicTransferEntropy(symbol, symbol, delta, numSymbols) == 0.0);
}

TEST_CASE("Test transfer entropy non_identical") {
    using namespace fast_entropy;
    int numSymbols = 4;
    int numLevels = 3;
    int maxExpectedNodes = 1000;

    SymbolSequence symbol1 = SymbolSequence(1000);
    SymbolSequence symbol2 = SymbolSequence(1000);
    for (int i=0; i < 1000; i++){
        int rand1 = rand() % numSymbols;
        int rand2 = rand() % numSymbols;
        symbol1[i] = rand1;
        symbol2[i] = rand2;
    }

    float transferEntropy = calculateSymbolicTransferEntropy(symbol1, symbol2, numLevels, numSymbols);
    REQUIRE(fabs(transferEntropy) < .05);
}
