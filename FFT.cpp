using namespace std;
#define _USE_MATH_DEFINES
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <complex>

const double TAU = 2. * M_PI;
const complex<double> IM_UNIT(0., 1.);
const complex<double> I_TAU(0., TAU);

static unsigned int num_digits (unsigned int array_len) {
    unsigned int digits_len = 0;
    for (int digit = array_len; digit > 0; digits_len++) {
        digit = digit >> 1;
    }
    --digits_len;
    cout << "== Array of length " << array_len;
    cout << " has " << digits_len << " address digits.\n";
    return digits_len;
}

static vector<complex<double>> bit_reverse_sort(vector<complex<double>>& in_array, unsigned int digits_len) {
    unsigned int array_len = in_array.size();

    vector<complex<double>> out_array { in_array };
    for (unsigned int idx = 0; idx < array_len; idx++) {
        unsigned int new_idx = 0;
        for (unsigned int digit = 0; digit < digits_len; digit++) {
            new_idx |= ((idx >> (digits_len - digit - 1)) & 1) << digit;
        }
        //cout << "item at index " << idx << " moved to index " << new_idx << "\n";
        out_array[new_idx] = in_array[idx];
    }
    return out_array;
}

int main() {
    // set up an array of complex numbers
    vector<complex<double>> array;
    unsigned int array_len{ 256 };               // <--- set desired input size here (order of 2)
    array.reserve(array_len);
    // get the length of the index int for array of this length
    const unsigned int digits_len{ num_digits(array_len) };

    // generate some input data
    float fl_arr_len = float(array_len - 1);
    for (unsigned int i = 0; i < array_len; i++) {
        array.emplace_back(complex<double>(cos(40. * TAU + i),0.));
        //array.emplace_back(complex<double>(float(i), float(fl_arr_len - i)));
    }

    /* print the input array to the terminal
    cout << "Initial array:\n";
    for (unsigned int i = 0; i < array.size(); i++) { 
        cout << "item " << i << ": " << array[i] << "\n"; 
    }
    cout << "Input array size: " << array.size() << ", and capacity: "
         << array.capacity() << "\n";*/

    // print the input data to a file
    ofstream input_file_out("fft_input.txt");

    for (complex<double> item : array) {
        input_file_out << item << "\n";
    }

    input_file_out.close();

    // Prepare for FFT with even/odd recursive sort
    vector<complex<double>> new_array = bit_reverse_sort(array, digits_len);
    cout << "Array from bit reverse size: " << new_array.size() << ", and capacity: "
         << new_array.capacity() << "\n";

    // radix-2 FFT
    for (unsigned int stage = 2; stage <= (array_len+1); stage *= 2) {
        // stage is 'm', the size of blocks, stage_size is 's'
        unsigned int stage_size = unsigned int(array_len / stage);
        cout << "m = " << stage << ", s = " << stage_size << "\n";

        complex<double> twiddle_scale_factor( exp(-I_TAU / double(stage) ));

        for (unsigned int block_id = 0; block_id < stage_size; block_id += stage) {
            // block_id is 'k' - inc by m to n-1
            cout << "k = " << block_id << "\n";

            complex<double> twiddle(1., 0.);

            for (unsigned int idx = 0; idx < (stage/2); idx++) {
                // idx is 'j' - inc to m/2
                // 2-pt butterfly (process pairs of data in blocks)
                unsigned int thisside = block_id + idx;                     // yk0 = j + k
                unsigned int thatside = thisside + stage / 2;               // yk1 = j + k + m//2
                cout << "ins: " << new_array[thisside] << ", " << new_array[thatside] << "\n";
                complex<double> t = twiddle * new_array[thatside];          // t = w * yk1
                
                new_array[thatside] = complex<double>(new_array[thisside] - t); // yk1 = u - t
                
                new_array[thisside] = complex<double>(new_array[thisside] + t); // yk0 = u + t

                cout << "outs: " << new_array[thisside] << ", " << new_array[thatside] << "\n";
                twiddle *= twiddle_scale_factor;
            }
        }
    }

    /* print the output array to the terminal
    cout << "Output array:\n";
    for (unsigned int i = 0; i < new_array.size(); i++) { 
        cout << "item " << i << ": " << new_array[i] << "\n";
    }
    cout << "Output array size: " << new_array.size() << ", and capacity: " 
         << new_array.capacity() << "\n";*/

    // print the output array to a file
    ofstream result_file_out("cpp_fft_out.txt");

    for (complex<double> item: new_array) {
        result_file_out << item << "\n";
    }

    result_file_out.close();

    return 0;
}
