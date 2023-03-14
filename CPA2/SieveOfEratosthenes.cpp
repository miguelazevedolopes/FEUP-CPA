#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>


void SieveOfEratosthenes(long long n, std::ofstream &outputFile)
{
    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long size = n - 1;
    bool *marks = new bool[size];

    std::cout << "Calculating SieveOfEratosthenes" << std::endl;

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
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
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

void SieveOfEratosthenesFastMarking(long long n, std::ofstream &outputFile)
{
    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long size = n - 1;
    bool *marks = new bool[size];

    std::cout << "Calculating SieveOfEratosthenesFastMarking" << std::endl;

    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    // Calculate primes
    for (long long k = 0; k < size; k++)
    {
        // Number is always 2 + index of array
        long long num = k + 2;

        // Go from num * num until n to check for multiples
        long long counter = 1;
        for (long long multNum = num * num; multNum <= n; counter++)
        {
            // Mark number if divisible by study number
            marks[multNum - 2] |= (multNum % num == 0);

            // Next number is always the base num * num plus the multiple of the base number
            multNum = num * num + num * counter;
        }
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
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

// TODO
void SieveOfEratosthenesFastMarkingReorganized(long long n, std::ofstream &outputFile)
{
    
    

    std::cout << "Calculating SieveOfEratosthenesFastMarking" << std::endl;


    // Start time
    auto t1 = std::chrono::high_resolution_clock::now();

    // Array in style -> [2, 3, 4, 5, 6, 7, 8, ..., n]
    // Size is n - 1
    long long limit = floor(sqrt(n)) + 1;

    std::vector<long long> prime;
    prime.reserve(limit);

    bool *marks = new bool[n-1];

    // Calculate primes
    for (long long k = 0; k < limit; k++)
    {
        // Number is always 2 + index of array
        long long num = k + 2;

        // Go from num * num until n to check for multiples
        long long counter = 1;
        for (long long multNum = num * num; multNum <= n; counter++)
        {
            // Mark number if divisible by study number
            marks[multNum - 2] |= (multNum % num == 0);

            // Next number is always the base num * num plus the multiple of the base number
            multNum = num * num + num * counter;
        }
    }

    // Output Primes
    for (long long k = 0; k < limit; k++)
    {
        // Unmarked are prime
        if (!marks[k])
        {
            // Number is always 2 + index of array
            prime.push_back(k+2);
        }
    }


    long long low = limit;
    long long high = 2*limit;


    while (low < n)
    {
        if (high >= n)
           high = n;


        bool mark[limit+1];
 

        for (long long i = 0; i < prime.size(); i++)
        {

            long long loLim = floor(low/prime[i]) * prime[i];
            if (loLim < low)
                loLim += prime[i];
 

            for (long long j=loLim; j<high; j+=prime[i])
                mark[j-low] = false;
        }

        for (int i = low; i<high; i++)
            if (mark[i - low] == true)
                marks[i] = true;

        low = low + limit;
        high = high + limit;
    }

    auto t2 = std::chrono::high_resolution_clock::now();

    // Output Execution time
    outputFile << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
               << std::endl;
    std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0 << "s" << std::endl
              << std::endl;

    // Output Primes
    for (long long k = 0; k < n-1; k++)
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
    long long n;

    std::cout << "N: ";
    std::cin >> n;  
    std::cout << std::endl;

    // std::ofstream outputFile1;
    // outputFile1.open("sieveoferatosthenes.txt");
    // SieveOfEratosthenes(n, outputFile1);
    // outputFile1.close();

    std::ofstream outputFile2;
    outputFile2.open("sieveoferatosthenesfastmarking.txt");
    SieveOfEratosthenesFastMarking(n, outputFile2);
    outputFile2.close();

    std::ofstream outputFile3;
    outputFile3.open("sieveoferatosthenesfastmarkingreorganized.txt");
    SieveOfEratosthenesFastMarkingReorganized(n, outputFile3);
    outputFile3.close();
}