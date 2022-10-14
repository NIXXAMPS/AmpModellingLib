// AmpModelling.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <stdio.h>
#include "AmpModelling.h"


int main()
{
    int i;
    float output[44100];

    BQFilter f1(FilterType::HighPass, 44100, 1e3f, 0.707f, 1, 1);

    for (i = 0; i < 44100; i++) {
        output[i] = f1.run(1);
        if(i > 44090)
            std::cout << output[i] << "\n";
    }

    

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
