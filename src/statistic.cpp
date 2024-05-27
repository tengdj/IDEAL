#include <iostream>
#include <fstream>
#include <set>

int main() {
    std::ifstream inputFile("output2.txt");
    std::ofstream outputFile("unique2.txt");

    if (!inputFile.is_open()) {
        std::cerr << "Error opening input file!" << std::endl;
        return 1;
    }

    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }

    std::set<unsigned int> uniqueData;
    unsigned int number;

    // Read data from input file and insert into set
    while (inputFile >> number) {
        uniqueData.insert(number);
    }

    inputFile.close();

    // Write unique data to output file
    for (const auto& data : uniqueData) {
        outputFile << data << std::endl;
    }

    outputFile.close();

    std::cout << "Data de-duplication completed successfully." << std::endl;
    return 0;
}
