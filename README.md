# Fast Transfer Entropy

Fast Transfer Entropy is a single-header C++ library for calculating transfer entropy and related measures using a custom tree based on the probability of symbol sequences. The library is designed to be efficient and easily parallelizable for multiple time series.

## Theory
Transfer entropy is a measure of the amount of information transferred from one time series to another. It quantifies the influence of one time series on another by measuring the reduction in uncertainty about the future of one time series given the past of another. Specifically, we calculate the symbolic transfer entropy $TE(X \to Y)$ from time series $X$ to $Y$ as follows:

$$TE(X \to Y) = H(Y_{t} | Y_t^{(\delta)}, X_{t}^{(\delta)}) - H(Y_{t+\delta} | Y_{t}^{(\delta)}) \\ = \sum_{\Omega} p(y, y^{(\delta)}, x^{(\delta)}) \log \left( \frac{p(y_{t} | y_t^{(\delta)}, x_t^{(\delta)})}{p(y_{t} | y_t^{(\delta)})} \right)$$

Where $`Y^{(\delta)}_t$`$ is the delay vector given by $`\{Y_{t-1},... ,Y_{t-\delta+1}\}`$.

## Features

- Efficient calculation of symbolic transfer entropy and KL transfer entropy.
- Probability tree-based implementation for symbol sequences.
- Single-header library for easy integration.
- easily parallelizable for multiple time series.

## Installation

Simply include the [fast_transfer_entropy.hpp](include/fast_transfer_entropy.hpp) header file in your project:

```cpp
#include "fast_transfer_entropy.hpp"
```

## Usage

The library provides a function for calculating symbolic transfer entropy. The following example demonstrates how to calculate symbolic transfer entropy between two symbol sequences:

```cpp
#include "fast_transfer_entropy.hpp"
int main() {
  fast_entropy::SymbolSequence X = {0, 1, 2, 3, 0, 1};
  fast_entropy::SymbolSequence Y = {1, 2, 3, 0, 1, 2};
  int delta = 2;
  int numSymbols = 4;

  double transferEntropy =
    fast_entropy::calculateSymbolicTransferEntropy(X, Y, delta, numSymbols);
  return 0;
}
```

## Documentation
The documentation is available via doxygen. To generate the documentation, run the following commands:
```bash
doxygen Doxyfile
```
The documentation will be generated in the `docs` directory. Open the `index.html` file in a web browser to view the documentation.

## Tests
Unit tests are provided in the `tests` directory. To run the tests, compile the test files and execute the resulting binaries:

```bash
cmake -Bbuild
cmake --build build
./build/tests
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
