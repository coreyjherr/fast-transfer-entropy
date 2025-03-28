#include <catch2/catch_test_macros.hpp>
#include "../include/fast_transfer_entropy.hpp"

TEST_CASE("Tree insertions make sense")
{
    using namespace fast_entropy;
    int numSymbols = 4;
    int numLevels = 3;
    int maxExpectedNodes = 1000;

    ProbabilityTree tree(numSymbols, numLevels, maxExpectedNodes);

    SymbolSequence symbol1 = {0, 1, 2};
    SymbolSequence symbol2 = {1, 2, 3};

    tree.insertSymbol(symbol1);
    tree.insertSymbol(symbol2);

    // Check if the symbols were inserted correctly
    tree.updateBaseCount();
    // Check the base count for the first symbol
    REQUIRE(tree.rootNodes[0]->count == 1);
    REQUIRE(tree.rootNodes[1]->count == 1);

    // Check the total count
    REQUIRE(tree.totalCount == 2);

    // Check the leaf list size
    REQUIRE(tree.leafNodes.size() == 2);

    // Check the daughter symbols to make sure they match
    REQUIRE(tree.leafNodes[0]->symbol == 2);

}


