#include "iostream"
#include "../include/data_struct.h"

void test_data_read_write()
{
    data_struct mystruct;

    mystruct.dump_to_csv("test.txt", {{"run_tag", "-"}, {"Temperature", "20"}, {"Voltage", "11"}, {"mirror_scan", "20"}, {"momentum_scan", "10"}},
                         {{"1234", {{"Temperature", "555"}, {"Voltage", "55"}}},
                          {"2345", {{"Temperature", "666"}, {"Voltage", "55"}}},
                          {"3456", {{"Temperature", "777"}, {"Voltage", "55"}}},
                          {"4567", {{"Temperature", "888"}, {"Voltage", "55"}}},
                          {"5678", {{"Temperature", "999"}, {"Voltage", "55"}}}});
}

void read_runbook(const char *filename = "../data/Run Book [2024][dRICH][testbeam] Run Book - RunBook.csv")
{

    std::ifstream file_date(filename);

    if (!file_date.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return; // Exit the program with an error
    }
    string line;
    int i = 0;
    getline(file_date, line);
    getline(file_date, line);
    // getline(file_date,line);
    // getline(file_date,line);

    std::istringstream iss_key(line);
    std::string field_key;
    std::vector<std::pair<std::string, std::string>> keys;

    // reading the second line for the key

    getline(file_date, line); // remove the 3rd line
    while (std::getline(iss_key, field_key, ','))
    {
        if (!field_key.empty())
            keys.push_back({field_key, "-"}); // second value is for default value will be updated later
    }
    int counter = 0; // just to limit for test study. removed later
    // reading rest of the data
    // std::string old_value = ""
    // std::vector<vector<string>>;

    // dump_to_csv(std::string output_file_name, std:: vector<pair<std::string, std::string>> data_structure, std::vector<pair<std::string, std::vector<std::pair<std::string, std::string>>>> data_dump = {});
    int m = 0;
    // for(auto key:keys){
    // cout << key.first << " "<< m++ <<endl;
    // }

    std::vector<std::pair<string, std::vector<std::pair<string, string>>>> dump_data;
    // std::vector <pair<string,string>> my
    while (std::getline(file_date, line))
    {

        std::istringstream iss(line);
        std::string field;
        int pos_counter = 0;
        // std::vector<string> line_val;
        // std::string
        std::string first_value = "norun";                   // first value of the data of dump_to_csv
        std::vector<std::pair<string, string>> second_value; // second value of the data of dump_to_csv
        while (std::getline(iss, field, ','))
        {

            // std::<pair> line_val;

            if (pos_counter > keys.size() - 1)
                break;

            // we have to make sure that if the key is run tag we do not re fill as each run number is the run tag.
            //  manually run tag entries are "calibrration , 15", "DCR scan" ,16 and Physics, 17;
            if (pos_counter != 15 || pos_counter != 16 || pos_counter == 17)
            {
                // cout << "run tag "<< keys[pos_counter].first <<endl;
                // continue;
                // }
                // else{
                if (!field.empty())
                    keys[pos_counter].second = field;
                if (field.empty())
                    field = keys[pos_counter - 1].second;

                second_value.push_back({keys[pos_counter].second, field});
            }

            if (pos_counter == 15 && !field.empty())
                first_value = field; // keys[pos_counter].first;
            else if (pos_counter == 16 && !field.empty())
                first_value = field; // keys[pos_counter].first;
            else if (pos_counter == 17 && !field.empty())
                first_value = field; // keys[17].first;
            // else continue;

            // line_val.push_back(field);

            // cout<< field << " ,";
            pos_counter++;

        } // with in line
        cout << first_value << "    " << counter << endl;
        // cout << line_val.size()<<endl;
        dump_data.push_back({first_value, second_value});
        // cout<<"\n";
        counter++;
        // if (counter>3) return;
    } // while rest of the data

    cout << dump_data[0].first << endl;

    data_struct mystruct2;
    mystruct2.dump_to_csv("result.txt", keys, dump_data);

} // end read runbook