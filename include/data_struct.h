#pragma once
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>

class data_struct
{

public:
    data_struct()
    {
    } // const
    ~data_struct()
    {

    } // dest
    void dump_to_csv(std::string output_file_name, std::vector<pair<std::string, std::string>> data_structure, std::vector<pair<std::string, std::vector<std::pair<std::string, std::string>>>> data_dump = {});

    void read_csv(const char *filename);

private:
    float xpos = 3.0;
};

void data_struct::dump_to_csv(std::string output_file_name, std::vector<pair<std::string, std::string>> data_structure, std::vector<pair<std::string, std::vector<std::pair<std::string, std::string>>>> data_dump = {})
{

    //  Open output stream
    std::ofstream output_file(output_file_name);

    //  Time stamp on creation
    std::time_t current_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    //  Vector of lines, start with creation timestamp
    std::vector<std::string> all_lines = {std::string("# file created on ") + std::string(std::ctime(&current_time))};

    //  Category line
    all_lines.push_back("# ");
    for (auto [category, default_value] : data_structure)
        all_lines[1] += std::string(category) + std::string(",  ");
    all_lines[1] = all_lines[1].substr(0, all_lines[1].size() - 3);
    all_lines[1] += std::string("\n");
    //   Data line
    auto current_pos = 1;
    for (auto [key, data_vect] : data_dump)
    {
        // cout<< key<< endl;
        all_lines.push_back(key + ",  ");
        current_pos++;
        // adding default value only if data vector is smaller than data structure size -1 (first is run tag)
        if (data_structure.size() - 1 > data_vect.size())
        {
            std::vector<int> nokey; // vector to store the position with no data
            // loop over data structure from 1 as zero is run tag
            for (int dat_i = 1; dat_i < data_structure.size(); dat_i++)
            {
                bool has_key = false; // bool to see if the data structure is present in data

                // loop over the data
                for (int dat_j = 0; dat_j < data_vect.size(); dat_j++)
                {

                    if (std::strcmp(data_vect[dat_j].first.c_str(), data_structure[dat_i].first.c_str()) == 0)
                        has_key = true; // if found the structure

                } // dat_j loop

                if (!has_key)
                {
                    nokey.push_back(dat_i);

                    data_vect.push_back({"", ""});
                }
            } // end data loop

            // loop to fill the default value of key from structure if not found in data
            for (auto key_pos : nokey)
            {
                data_vect[key_pos - 1].first = data_structure[key_pos].first;   // key
                data_vect[key_pos - 1].second = data_structure[key_pos].second; // default value

            } // end data structure position

        } // END if data_structure is greater than data vector provided
        for (auto [category, value] : data_vect)
        {

            all_lines[current_pos] += std::string(value) + std::string(",  ");
        }

        all_lines[current_pos] = all_lines[current_pos].substr(0, all_lines[current_pos].size() - 3);

        all_lines[current_pos] += std::string("\n");
    }
    //  Populate lines
    for (auto current_line : all_lines)
        output_file << current_line;
    output_file.close();
} // end dumb_to_csv

void read_csv(const char *filename)
{
    std::cout << filename << endl;
}