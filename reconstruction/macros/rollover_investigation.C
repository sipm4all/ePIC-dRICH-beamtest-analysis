#pragma once

class decoded_data
{
public:
    // Getters
    int get_device() const { return device; }
    int get_fifo() const { return fifo; }
    int get_type() const { return type; }
    int get_counter() const { return counter; }
    int get_column() const { return column; }
    int get_pixel() const { return pixel; }
    int get_tdc() const { return tdc; }
    int get_rollover() const { return rollover; }
    int get_coarse() const { return coarse; }
    int get_fine() const { return fine; }
    TTree *get_data_tree() const
    {
        cout << "Sure, here is your tree with " << data_tree->GetEntries() << " entries!" << endl;
        return data_tree;
    }
    TGraph *get_rollover_counter() const { return rollover_counter; }
    TH1F *get_h_counters() const { return h_counters; }

    // Setters
    void set_device(int value) { device = value; }
    void set_fifo(int value) { fifo = value; }
    void set_type(int value) { type = value; }
    void set_counter(int value) { counter = value; }
    void set_column(int value) { column = value; }
    void set_pixel(int value) { pixel = value; }
    void set_tdc(int value) { tdc = value; }
    void set_rollover(int value) { rollover = value; }
    void set_coarse(int value) { coarse = value; }
    void set_fine(int value) { fine = value; }
    void set_data_tree(TTree *value) { data_tree = value; }
    void set_rollover_counter(TGraph *value) { rollover_counter = value; }
    void set_h_counters(TH1F *value) { h_counters = value; }
    void set_variables_to_tree(TTree *tree);

    //  File I/O
    void load_data_file(std::string target_file);

protected:
private:
    int device = 0;
    int fifo = 0;
    int type = 0;
    int counter = 0;
    int column = 0;
    int pixel = 0;
    int tdc = 0;
    int rollover = 0;
    int coarse = 0;
    int fine = 0;
    TTree *data_tree = nullptr;
    TGraph *rollover_counter = nullptr;
    TH1F *h_counters = nullptr;
};
void decoded_data::set_variables_to_tree(TTree *tree)
{
    if (!tree)
        throw std::runtime_error("[ERROR] TTree pointer is null.");

    // Map TTree branches to class variables
    tree->SetBranchAddress("device", &device);
    tree->SetBranchAddress("fifo", &fifo);
    tree->SetBranchAddress("type", &type);
    tree->SetBranchAddress("counter", &counter);
    tree->SetBranchAddress("column", &column);
    tree->SetBranchAddress("pixel", &pixel);
    tree->SetBranchAddress("tdc", &tdc);
    tree->SetBranchAddress("rollover", &rollover);
    tree->SetBranchAddress("coarse", &coarse);
    tree->SetBranchAddress("fine", &fine);
}

void decoded_data::load_data_file(std::string target_file)
{
    // Open the file
    std::unique_ptr<TFile> input_file(TFile::Open(target_file.c_str(), "READ"));
    if (!input_file || input_file->IsZombie())
        throw std::runtime_error("[ERROR] Failed to open file: " + target_file);

    // Get components
    auto current_data_tree = dynamic_cast<TTree *>(input_file->Get("alcor"));
    if (!current_data_tree)
        throw std::runtime_error("[ERROR] TTree 'alcor' not found in file: " + target_file);

    // Clone the tree
    data_tree = current_data_tree->CloneTree();
    if (!data_tree)
        throw std::runtime_error("[ERROR] Failed to clone the TTree.");

    rollover_counter = dynamic_cast<TGraph *>(input_file->Get("gRollover"));
    if (!rollover_counter)
        throw std::runtime_error("[ERROR] TGraph 'gRollover' not found in file: " + target_file);

    h_counters = dynamic_cast<TH1F *>(input_file->Get("hCounters"));
    if (!h_counters)
        throw std::runtime_error("[ERROR] TH1F 'hCounters' not found in file: " + target_file);

    // Set tree to data
    try
    {
        set_variables_to_tree(data_tree);
    }
    catch (const std::exception &e)
    {
        std::cerr << "[ERROR] Exception occurred while setting TTree branches: " << e.what() << std::endl;
    }
}

void rollover_investigation()
{
    // Target file
    std::string target_file = "/Users/nrubini/Analysis/ePIC/_production_repositories/sipm4eic-testbeam2023-analysis/Data/20231010-084623/kc705-193/decoded/alcdaq.fifo_0.root";

    // Load data
    decoded_data data_container;
    data_container.load_data_file(target_file);
    auto current_tree = data_container.get_data_tree();

    // Loop over data
    for (auto i_ter = 0; i_ter < current_tree->GetEntries(); i_ter++)
    {
        //  Get entries
        current_tree->GetEntry(i_ter);

        // Example: Print the device value
        std::cout << "Entry " << i_ter << ", Device: " << data_container.get_device() << std::endl;
    }
}