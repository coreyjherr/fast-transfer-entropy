#include <catch2/catch_test_macros.hpp>
#include "../include/fast_transfer_entropy.hpp"

TEST_CASE("Test Entropy"){
    using namespace fast_entropy;
    int numSymbols = 4;
    int numLevels = 3;
    int maxExpectedNodes = 1000;

    ProbabilityTree tree(numSymbols, numLevels, maxExpectedNodes);
    for (int i=0; i < 100; i++){
        SymbolSequence symbol1 = {0, 0, 0};
        tree.insertSymbol(symbol1);
    }

    tree.updateBaseCount();
    REQUIRE(tree.rootNodes[0]->count == 100);

    // checkt the entropy is 0
    double entropy = tree.calculateEntropy();
    REQUIRE(entropy == 0);
}

TEST_CASE("Test Entropy of random symbols"){
    using namespace fast_entropy;
    for (int numSymbols = 3; numSymbols < 10; numSymbols++) {
        int numLevels = 3;
        int maxExpectedNodes = 10000;

        ProbabilityTree tree(numSymbols, numLevels, maxExpectedNodes);

        // Insert random symbols
        for (int i=0; i < 10000; i++){
            SymbolSequence symbol1 = {rand() % numSymbols, rand() % numSymbols, rand() % numSymbols};
            tree.insertSymbol(symbol1);
        }

        // Calculate entropy
        double entropy = tree.calculateEntropy();

        // Check if entropy is near the expected value of log2(numSymbols)
        REQUIRE(entropy-log2(numSymbols) < 0.05);
    }

}


