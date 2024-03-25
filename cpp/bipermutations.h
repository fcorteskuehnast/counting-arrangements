#pragma once
#include <vector>
#include <iostream>

const std::vector<std::string> knownValues = {
    "1",
    "1",
    "1",
    "2",
    "8",
    "62",
    "908",
    "24698",
    "1232944",
    "112018190",
    "18410581880",
    "5449192389984",
    "2894710651370536",
    "2752596959306389652",
    "4675651520558571537540",
    "14163808995580022218786390",
    "76413073725772593230461936736"
};

const std::vector<size_t> hexBip2 = {
    0, 1, 2, 3, 4, 5, 6,
    7, 8, 3, 9, 10, 11, 12,
    13, 4, 14, 15, 1,
    12, 5, 16, 17, 9, 18, 0, 6,
    11, 15, 2, 17, 13,
    7, 18, 10, 14, 16, 8 };

const std::vector<size_t> hexBip2x3x3 = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 8, 14, 15, 16,
    9, 17, 18, 6, 19, 10,
    3, 20, 14, 21, 5, 11, 22, 16, 7,
    0, 19, 4, 12, 17, 1,
    20, 15, 22, 2, 13, 18, 21 };

// midpoint is 3-wise intersection, contains 7 6-wise intersections
const std::vector<size_t> hexBip_m3_x7 = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 11, 9, 12, 13, 14, 15, 16,
    5, 17, 18, 19, 20, 21, 22, 23, 6,
    15, 26, 19, 10, 0, 24, 7, 21, 14, 4,
    25, 18, 1, 8, 16, 22, 13, 24,
    25, 26, 17, 11, 23, 2, 3, 12, 20 };

const std::vector<size_t> hexBip3 = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 8, 14, 15, 16, 17, 18,
    9, 19, 20, 21, 6, 22, 23, 3, 24,
    10, 25, 26, 17, 27, 5, 14, 28, 21, 11, 0,
    24, 7, 16, 22, 4, 26, 20, 12, 1,
    18, 28, 23, 15, 27, 2, 13, 19, 25 };

// hexagonal bipermutations in the region with 5 slopes
const std::vector<size_t> hexBip3_5_slopes = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 8, 13, 14, 15, 9,
    16, 17, 6, 18, 3,
    10, 19, 20, 15, 21, 5, 13, 22, 17, 11, 0,
    7, 14, 4, 20, 12, 1,
    22, 18, 21, 2, 16, 19 };

// bipermutations of kxl rectangles in R_4
const std::vector<size_t> rectBip2x2 = { 0, 2, 1, 4, 3, 5, 6, 1, 3, 0, 5, 2, 4, 6 };
const std::vector<size_t> rectBip2x3 = { 0, 2, 1, 4, 3, 5, 7, 6, 1, 8, 0, 3, 6, 8, 5, 2, 4, 7 };
const std::vector<size_t> rectBip3x4 = { 0, 3, 1, 10, 2, 4, 6, 10, 7, 5, 8, 11, 9, 2, 12, 1, 5, 0, 9, 4, 8, 12, 7, 3, 6, 11 };

const std::vector<size_t> rectBip4x4 = { 0, 4, 1, 12, 2, 5, 3, 8, 6, 9, 12, 10, 7, 11, 13, 3, 7, 2, 14, 1, 6, 0, 11, 5, 10, 14, 9, 4, 8, 13 };
const std::vector<size_t> rectBip4x5 = { 0, 4, 1, 13, 2, 5, 3, 8, 6, 9, 13, 10, 7, 11, 14, 12, 3, 15, 2, 7, 1, 16, 0, 6, 12, 16, 11, 5, 10, 15, 9, 4, 8, 14 };
const std::vector<size_t> rectBip4x7 = { 0, 4, 1, 16, 2, 5, 3, 9, 6, 10, 16, 11, 7, 12, 17, 13, 8, 14, 18, 15, 3, 19, 2, 8, 1, 20, 0, 7, 15, 20, 14, 6, 13, 19, 12, 5, 11, 18, 10, 4, 9, 17 };
const std::vector<size_t> rectBip4x9 = { 0, 4, 1, 19, 2, 5, 3, 10, 6, 11, 19, 12, 7, 13, 20, 14, 8, 15, 21, 16, 9, 17, 22, 18, 3, 23, 2, 9, 1, 24, 0, 8, 18, 24, 17, 7, 16, 23, 15, 6, 14, 22, 13, 5, 12, 21, 11, 4, 10, 20 };
const std::vector<size_t> rectBip5x6 = { 0, 5, 1, 16, 2, 6, 20, 3, 4, 7, 10, 20, 11, 8, 12, 16, 13, 9, 14, 17, 15, 4, 18, 3, 9, 2, 19, 1, 8, 0, 15, 7, 14, 19, 13, 6, 12, 18, 11, 5, 10, 17 };

const std::vector<size_t> rectBip6x6 = { 0, 6, 1, 19, 2, 7, 3, 18, 4, 8, 5,
                                        12, 9, 13, 18, 14, 10, 15, 19, 16, 11, 17, 20,
                                        5, 11, 4, 21, 3, 10, 2, 22, 1, 9, 0,
                                        17, 8, 16, 22, 15, 7, 14, 21, 13, 6, 12, 20 };

const std::vector<size_t> rectBip6x7 = { 0, 6, 1, 20, 2, 7, 3, 19, 4, 8, 5,
                                        12, 9, 13, 19, 14, 10, 15, 20, 16, 11, 17, 21, 18,
                                        5, 22, 4, 11, 3, 23, 2, 10, 1, 24, 0,
                                        9, 18, 24, 17, 8, 16, 23, 15, 7, 14, 22, 13, 6, 12, 21 };

const std::vector<size_t> rectBip7x8 = { 0, 7, 1, 24, 2, 8, 3, 23, 4, 9, 5, 22, 6,
                                        10, 14, 22, 15, 11, 16, 23, 17, 12, 18, 24, 19, 13, 20, 25, 21,
                                        6, 26, 5, 13, 4, 27, 3, 12, 2, 28, 1, 11, 0,
                                        21, 10, 20, 28, 19, 9, 18, 27, 17, 8, 16, 26, 15, 7, 14, 25 };

const std::vector<size_t> rectBip8x8 = { 0, 8, 1, 26, 2, 9, 3, 25, 4, 10, 5, 24, 6, 11, 7,
                                        16, 12, 17, 24, 18, 13, 19, 25, 20, 14, 21, 26, 22, 15, 23, 27,
                                        7, 15, 6, 28, 5, 14, 4, 29, 3, 13, 2, 30, 1, 12, 0,
                                        23, 11, 22, 30, 21, 10, 20, 29, 19, 9, 18, 28, 17, 8, 16, 27 };

const std::vector<size_t> rectBip8x9 = { 0, 8, 1, 27, 2, 9, 3, 26, 4, 10, 5, 25, 6, 11, 7,
                                        16, 12, 17, 25, 18, 13, 19, 26, 20, 14, 21, 27, 22, 15, 23, 28, 24,
                                        7, 29, 6, 15, 5, 30, 4, 14, 3, 31, 2, 13, 1, 32, 0,
                                        12, 24, 32, 23, 11, 22, 31, 21, 10, 20, 30, 19, 9, 18, 29, 17, 8, 16, 28 };

std::vector<size_t> get_square_bip(int region, int size)
{
    /**
     * calculates the boundary bipermutation of a square that is aligned with the axes
     */
    std::vector<size_t> bip;
    if (region == 3)
    {
        for (int i = 0; i < size; i++)
        {
            bip.push_back(i);
            bip.push_back(4 * size - 2 - i);
        }
        bip.pop_back();

        for (int i = 0; i < size; i++)
        {
            bip.push_back(3 * size - 1 - i);
            bip.push_back(size + i);
        }

        for (int i = 0; i < size; i++)
        {
            bip.push_back(size - 1 - i);
            bip.push_back(2 * size + i);
        }
        bip.pop_back();

        for (int i = 0; i < size; i++)
        {
            bip.push_back(3 * size - 1 + i);
            bip.push_back(2 * size - 1 - i);
        }

        return bip;
    }
    else if (region == 4)
    {
        for (int i = 0; i < size; i++)
        {
            bip.push_back(5 * size - 2 - i);
            bip.push_back(i);
            bip.push_back(4 * size - 2 - i);
        }
        bip.pop_back();

        for (int i = 0; i < size; i++)
        {
            bip.push_back(3 * size - 1 - i);
            bip.push_back(size + i);
            bip.push_back(4 * size - 1 + i);
        }
        // bip.pop_back();

        for (int i = 0; i < size; i++)
        {
            // bip.push_back(5 * size - 2 + i);
            bip.push_back(size - 1 - i);
            bip.push_back(5 * size - 2 + 1 + i); //
            bip.push_back(2 * size + i);
        }
        // bip.pop_back();

        for (int i = 0; i < size; i++)
        {
            // bip.push_back(3 * size - 1 + i);
            // bip.push_back(2 * size - 1 - i);
            // bip.push_back(6 * size - 3 - i);
            bip.push_back(6 * size - 2 - i);
            bip.push_back(2 * size - 1 - i);
            bip.push_back(3 * size - 1 + 1 + i);
        }
        bip.pop_back();

        return bip;
    }

    std::cout << "region is not 3 or 4";
    throw std::invalid_argument("region is not 3 or 4");
}
