#include <iostream>
#include <fstream>
#include <chrono>

void SieveOfEratosthenes(long long n, std::ofstream &outputFile)
{
    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long size = n - 1;
    bool *marks = new bool[size];

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    // Calculate primes
    for (long long k = 0; k < size; k++)
    {
        // Number is always 2 + index of array
        long long num = k + 2;

        // Go from num * num until n to check for multiples
        for (long long j = num * num - 2; j < size; j++)
        {
            // Number is always 2 + index of array
            long long multNum = j + 2;

            // Mark number if divisible by study number
            marks[j] |= (multNum % num == 0);
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
               << std::endl;

    // Output Primes
    for (long long k = 0; k < size; k++)
    {
        // Unmarked are prime
        if (!marks[k])
        {
            // Number is always 2 + index of array
            outputFile << k + 2 << std::endl;
        }
    }

    // Clear memory :D
    delete marks;
}

int main(int argc, char *argv[])
{
    std::ofstream outputFile;
    outputFile.open("sieveoferatosthenes.txt");

    SieveOfEratosthenes(1000000, outputFile);
}