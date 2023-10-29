#include <iostream>
#include <fstream>

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const uint32_t un32 = 1;  // we only need the first bit here -->> will be used to check the lowest bit of other unsigned integer (with arbitrarily larger number of bits)

//const unsigned int n_max = 128;  // for bitset

/******************************************************************************/
/*******************   Convert Integer to Binary string   *********************/
/******************************************************************************/
std::string int_to_bstring(__int128_t bool_nb, unsigned int r)
{
    std::string s;
    do
    {
        s.push_back( ((bool_nb & un32)?'1':'0') );
    } while(bool_nb >>= 1);

    reverse(s.begin(), s.end());
    s = (std::string(r - s.length(), '0')).append(s);

    return s;
}

// Using bitset: Very bad performance (~ double the time)
/*std::string int_to_bstring_Bitset(__int128_t bool_nb)     // Very bad performance (~ double the time)
{
    std::bitset<n_max> hi{ static_cast<unsigned long long>(bool_nb >> 64) },
            lo{ static_cast<unsigned long long>(bool_nb) },
            bits{ (hi << 64) | lo };
    return bits.to_string();
}*/

/******************************************************************************/
/****************   Count number of set bits of an integer  *******************/
/******************************************************************************/
unsigned int bitset_count(__int128_t bool_nb)      // Using trick; by far the fastest option in all cases (faster than using bitset or checking bits one by one)
{
    unsigned int count;
    for (count=0; bool_nb; count++)
        bool_nb &= bool_nb - 1;
    return count;
}

/******************************************************************************/
/**********************   Position of the lowest bit  *************************/
/******************************************************************************/
// returns an integer with only one single bit at the position of the lowest bit of "a"
__int128_t Lowest_Bit(__int128_t a)
{
    return ((a - 1) ^ a) & a;
}

/******************************************************************************/
/****************   Print the digit location of each bits   *******************/
/******************************************************************************/
const unsigned int un128 = 1;

void int_to_digits(__int128_t bool_nb, unsigned int r)
{
    unsigned int digit = r;

    std::cout << "\t"; // << std::endl;
    while(bool_nb)
    {
        if(bool_nb & un128) {   std::cout << digit << "\t"; }
        bool_nb >>= 1;
        digit--;
    }
    std::cout << std::endl;
}


void int_to_digits_file(__int128_t bool_nb, unsigned int r, std::fstream &file)
{
    unsigned int digit = r;

    file << "\t"; // << std::endl;
    while(bool_nb)
    {
        if(bool_nb & un128) {   file << digit << "\t"; }
        bool_nb >>= 1;
        digit--;
    }
    file << std::endl;
}

/******************************************************************************/
/********************************   N CHOOSE K   ******************************/
/******************************************************************************/

unsigned int Choose(unsigned int n, unsigned int k)
{
    double res = 1;
    for (unsigned int i = 1; i <= k; ++i)
        res = res * (n - k + i) / ((double) i);
    return (unsigned int)(res + 0.01);
}
