#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <limits>
#include <deque>
#include <string>
#include <sstream>
#include <dirent.h>
#include <sys/stat.h>
#include <filesystem>
#include <chrono>
#include <unistd.h>
using namespace std;

const uint64_t N = 32;
const size_t size_alphavite = 256;
const uint64_t pow_of_N = (1ull << N);
const uint64_t midle_of_interval = (pow_of_N >> 1);
const uint64_t bin_pow_of_N_minus_2 = (pow_of_N >> 2);
size_t fsz=0;

bool _1 = false;
bool _9 = false;
bool _d = false;
bool _t = false;
bool _r = false;
bool _l = false;
bool _k = false;
bool _c = false;

struct coded_message {
    vector<bool> codded_string;
    size_t size;
};

uint64_t to_unsigned_shar(const char ch) {
    return ch + 128;
}

char to_signed_shar(const uint64_t ind) {
    return ind - 128;
}

void ar_send_bit(const bool bit, int64_t& waiting_counter,vector<bool>& codded_string) {
    codded_string.push_back(bit);
    while (waiting_counter > 0) {
        if (bit == 0) {
            codded_string.push_back(1);
        } else {
            codded_string.push_back(0);
        }
        --waiting_counter;
    }
}

pair<uint64_t, uint64_t> update_current_interval(const uint64_t l, const uint64_t h, const uint64_t ind, const uint64_t cumulitive_ind, const vector<uint64_t>& cum_freq) {
    return std::make_pair((l + (h - l + 1) * cum_freq[ind] / cum_freq[cumulitive_ind]), (l + (h - l + 1) * cum_freq[ind + 1] / cum_freq[cumulitive_ind] - 1));
}

coded_message compress(const string& text) {
    const size_t size = text.size();

    vector<uint64_t> cum_freq(size_alphavite + 1, 0);
    vector<uint64_t> freq(size_alphavite, 1);
    vector<bool> codded_string;
    cum_freq[0] = 0;
    for (size_t i = 0; i < size_alphavite; ++i) //���������� ������, ���� ������� ���� ������� cum_freq
        cum_freq[i + 1] = cum_freq[i] + freq[i];

    uint64_t l = 0, h = pow_of_N - 1;//2^N-1
    int64_t waiting_counter = 0;
    for (size_t i = 0; i < size; ++i) {
        size_t ind = to_unsigned_shar(text[i]);
        while (1) {
            if (h < midle_of_interval) {//���������� ������� ����� � ��������� 0(<10000000000)
                ar_send_bit(0, waiting_counter, codded_string);
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
            } else if (l >= midle_of_interval) {//���������� ������� ����� � ��������� 1
                ar_send_bit(1, waiting_counter, codded_string);
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
            }
            else if ((h - l) < cum_freq[size_alphavite]) {//�������� ������� ��������
                ++waiting_counter;
                l = (l - bin_pow_of_N_minus_2) << 1;
                h = ((h - bin_pow_of_N_minus_2) << 1) + 1;
            } else { break; }

        }
        pair<uint64_t, uint64_t> new_interval = update_current_interval(l, h, ind, size_alphavite, cum_freq);
        l = new_interval.first;
        h = new_interval.second;
        //update cum_freq � freq
        ++freq[ind];
        for (size_t j = ind; j < size_alphavite; ++j)
            ++cum_freq[j + 1];
        //����� ���������� �� ��������� - ������ �� ������� �  EOF - ���������� ���������� �������.
        if (i == size - 1) {
            for (size_t j = 0; j < N; ++j) {
                uint64_t cur_left_oldest_bit = l & midle_of_interval;
                ar_send_bit((cur_left_oldest_bit) >> (N - 1), waiting_counter, codded_string);
                l = l << 1;
            }
        }

    }
    return { codded_string, size };
}

string decompress(const coded_message& coded_message) {

    vector<uint64_t> cum_freq(size_alphavite + 1, 0);
    vector<uint64_t> freq(size_alphavite, 1);

    cum_freq[0] = 0;
    for (size_t i = 0; i < size_alphavite; ++i) //���������� ������, ���� ������� ���� ������� cum_freq
        cum_freq[i + 1] = cum_freq[i] + freq[i];

    string decoded_text;

    uint64_t value = coded_message.codded_string[0];
    size_t i = 1;
    for (; i < N; ++i) {
        value = (value << 1) + coded_message.codded_string[i];
    }
    uint64_t l = 0, h = pow_of_N - 1;//2^N-1

    for (size_t cur_sz = 0; cur_sz < coded_message.size; ++cur_sz) {
        //Decode 1 symbol
        for (size_t ind = 0; ind < size_alphavite; ++ind) {

            pair<uint64_t, uint64_t> new_interval = update_current_interval(l, h, ind, size_alphavite, cum_freq);
            uint64_t l_of_symb = new_interval.first;
            uint64_t h_of_symb = new_interval.second;

            if (l_of_symb <= value && value <= h_of_symb) {
                decoded_text.push_back(to_signed_shar(ind));
                l = l_of_symb;
                h = h_of_symb;
                ++freq[ind];
                for (size_t j = ind; j < size_alphavite; ++j)
                    ++cum_freq[j + 1];
                break;
            }

        }
        //
        while (1) {
            if (h < midle_of_interval) {//���������� ������� ����� � ��������� 0(<10000000000)
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
                value = ((value << 1) % pow_of_N) + coded_message.codded_string[i++];
            } else if (l >= midle_of_interval) {//���������� ������� ����� � ��������� 1
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
                value = ((value << 1) % pow_of_N) + coded_message.codded_string[i++];
            }
            else if ((h - l) < cum_freq[size_alphavite]) {//�������� ������� ��������

                l = (l - bin_pow_of_N_minus_2) << 1;
                h = ((h - bin_pow_of_N_minus_2) << 1) + 1;
                value = ((value - bin_pow_of_N_minus_2) << 1) + coded_message.codded_string[i++];
            } else { break; }

        }

    }

    return decoded_text;
}

void ar_send_bit(const bool bit, int64_t& waiting_counter, deque<bool>& codded_string) {
    codded_string.push_back(bit);
    while (waiting_counter > 0) {
        if (!bit) {
            codded_string.push_back(1);
        }
        else {
            codded_string.push_back(0);
        }
        --waiting_counter;
    }
}

void arif_compress(std::istream& fin, std::ostream& fout) {
    std::ofstream t_fout("temp_output.bin", ios::binary);
    if (!t_fout) {
        cerr << "Error creating temporary file." << endl;
        return;
    }
    vector<uint64_t> cum_freq(size_alphavite + 1, 0);
    vector<uint64_t> freq(size_alphavite, 1);
    deque<bool> codded_string;
    uint64_t l = 0, h = pow_of_N - 1;//2^N-1
    cum_freq[0] = 0;
    for (size_t i = 0; i < size_alphavite; ++i) //���������� ������, ���� ������� ���� ������� cum_freq
        cum_freq[i + 1] = cum_freq[i] + freq[i];

    int64_t waiting_counter = 0;
    char ch;
    uint64_t letters = 0;
    while (fin.get(ch)) {

        ++letters;

        size_t ind = to_unsigned_shar(ch);
        while (1) {
            if (h < midle_of_interval) {//���������� ������� ����� � ��������� 0(<10000000000)
                ar_send_bit(0, waiting_counter, codded_string);
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
                while (codded_string.size() >= 8) {
                    char c = 0;
                    for (int i = 0; i < 8; ++i) {
                        c |= (codded_string.front()) << (7 - i);
                        codded_string.pop_front();
                    }
                    t_fout.put(c);
                }

            } else if (l >= midle_of_interval) {//���������� ������� ����� � ��������� 1
                ar_send_bit(1, waiting_counter, codded_string);
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
                while (codded_string.size() >= 8) {
                    char c = 0;
                    for (int i = 0; i < 8; ++i) {
                        c |= (codded_string.front()) << (7 - i);
                        codded_string.pop_front();
                    }
                    t_fout.put(c);
                }
            }
            else if ((h - l) < cum_freq[size_alphavite]) {//�������� ������� ��������
                ++waiting_counter;
                l = (l - bin_pow_of_N_minus_2) << 1;
                h = ((h - bin_pow_of_N_minus_2) << 1) + 1;
            } else { break; }

        }
        pair<uint64_t, uint64_t> new_interval = update_current_interval(l, h, ind, size_alphavite, cum_freq);
        l = new_interval.first;
        h = new_interval.second;
        //update cum_freq � freq
        ++freq[ind];
        for (size_t j = ind; j < size_alphavite; ++j)
            ++cum_freq[j + 1];

    }

    for (size_t j = 0; j < N; ++j) {
        uint64_t cur_h_oldest_bit = h & midle_of_interval;
        ar_send_bit((cur_h_oldest_bit) >> (N - 1), waiting_counter, codded_string);
        h = h << 1;
    }

    while (codded_string.size() % 8 != 0) {
        codded_string.push_back(0);
    }
    while (codded_string.size() != 0) {
        char c = 0;
        for (int i = 0; i < 8; ++i) {
            c |= (codded_string.front()) << (7 - i);
            codded_string.pop_front();
        }
        t_fout.put(c);
    }
    t_fout.close();
    std::ifstream ti_fout("temp_output.bin", ios::binary);
    fout.write(reinterpret_cast<const char*>(&fsz), sizeof(fsz));
    fout.write(reinterpret_cast<const char*>(&letters), sizeof(letters));
    char c;
    while (ti_fout.get(c)) {
        fout.put(c);
    }

    ti_fout.close();

    if (remove("temp_output.bin") != 0) {
        cerr << "Error deleting temporary file." << endl;
    }
}

bool getBit(char ch, int bitPosition) {
    return (ch & (1 << bitPosition)) != 0;
}

void arif_decompress(std::istream& fin, std::ostream& fout) {
    vector<uint64_t> cum_freq(size_alphavite + 1, 0);
    vector<uint64_t> freq(size_alphavite, 1);

    cum_freq[0] = 0;
    for (size_t i = 0; i < size_alphavite; ++i) //���������� ������, ���� ������� ���� ������� cum_freq
        cum_freq[i + 1] = cum_freq[i] + freq[i];
    uint64_t value = 0;
    size_t i = 0;
    char ch, dec;
    uint64_t l = 0, h = pow_of_N - 1;//2^N-1
    uint64_t letters = 0;
    fin.read(reinterpret_cast<char*>(&fsz), sizeof(fsz));
    fin.read(reinterpret_cast<char*>(&letters), sizeof(letters));
    for (; i < N; ++i) {
        if (i % 8 == 0) {
            fin.get(ch);

        }
        bool biting = getBit(ch, (7 - i % 8));
        value = (value << 1) + biting;
    }
    for (size_t cur_sz = 0; cur_sz < letters; ++cur_sz) {

        for (size_t ind = 0; ind < size_alphavite; ++ind) {

            pair<uint64_t, uint64_t> new_interval = update_current_interval(l, h, ind, size_alphavite, cum_freq);
            uint64_t l_of_symb = new_interval.first;
            uint64_t h_of_symb = new_interval.second;

            if (l_of_symb <= value && value <= h_of_symb) {
                dec = to_signed_shar(ind);
                fout.put(dec);
                l = l_of_symb;
                h = h_of_symb;
                ++freq[ind];
                for (size_t j = ind; j < size_alphavite; ++j)
                    ++cum_freq[j + 1];
                break;
            }

        }

        //
        while (1) {
            if (h < midle_of_interval) {//���������� ������� ����� � ��������� 0(<10000000000)
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
                if (i % 8 == 0) {
                    fin.get(ch);

                }
                bool biting = getBit(ch, (7 - i % 8));
                value = ((value << 1) % pow_of_N) + biting;
                ++i;
            } else if (l >= midle_of_interval) {//���������� ������� ����� � ��������� 1
                l = (l << 1) % pow_of_N;//����� �� 1
                h = (h << 1) % pow_of_N + 1;
                if (i % 8 == 0) {
                    fin.get(ch);

                }
                bool biting = getBit(ch, (7 - i % 8));
                value = ((value << 1) % pow_of_N) + biting;
                ++i;
            }
            else if ((h - l) < cum_freq[size_alphavite]) {//�������� ������� ��������

                l = (l - bin_pow_of_N_minus_2) << 1;
                h = ((h - bin_pow_of_N_minus_2) << 1) + 1;
                if (i % 8 == 0) {
                    fin.get(ch);

                }
                bool biting = getBit(ch, (7 - i % 8));
                value = ((value - bin_pow_of_N_minus_2) << 1) + biting;
                ++i;
            } else { break; }

        }

    }

}

string read_File(const string& file_path) {

    ifstream file(file_path, ios::binary);

    if (!file) {
        cerr << "�� ������� ������� ����: " << file_path << endl;
        cout << "�� ������� ������� ����: " << file_path << endl;
        return "";
    }

    file.seekg(0, ios::end);
    streampos file_size = file.tellg();
    file.seekg(0, ios::beg);

    string s(file_size, 'a');

    if (!file.read(&s[0], file_size)) {
        cerr << "�� ������� ��������� ����: " << file_path << endl;
        cout << "�� ������� ��������� ����: " << file_path << endl;
        return "";
    }

    file.close();

    return s;
}

void writeStringToBinaryFile(const std::string& filename, const std::string& data) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error opening file for writing." << std::endl;
        return;
    }

    // Write the string data to the file
    outfile.write(data.c_str(), data.size());

    // Close the file
    outfile.close();
}

//LZ77
struct Node {
    int64_t offset;
    int64_t length;
    char next;
};

Node findMatching(string& s, int64_t pos, int64_t window_sz)
{

    Node ans = { 0, 0, '\0' };

    if (_9) { window_sz = pos; }
    for (int64_t i = max((int64_t)1, pos - window_sz); i <= pos; ++i) {
        int64_t j = 0;
        while (pos + j < s.size() && s[pos + j] == s[pos - i + j])
            ++j;

        if (j > ans.length) {
            ans.length = j;
            ans.offset = i;
        }
    }

    if (ans.length + pos < s.length()) {
        ans.next = s[pos + ans.length];
    }

    return ans;
}

void encodeLZ77(string& s, vector<Node>& ans)
{

    size_t pos = 0;
    int64_t window_sz = 4096;
    if(_1) window_sz=1024;
    Node temp;
    while (pos < s.length()) {

        temp = findMatching(s, pos, window_sz);

        pos += temp.length + 1;

        ans.push_back(temp);
    }
}

string decodeLZ77(vector<Node>& crypted)
{

    string ans = "";
    for (size_t i = 0; i < crypted.size(); ++i) {
        if (crypted[i].length > 0) {
            int64_t start = ans.size() - crypted[i].offset;
            for (int64_t j = 0; j < crypted[i].length; ++j) {

                ans += ans[start + j];
            }
        }

        ans += crypted[i].next;
    }
    return ans;
}

pair<string, string> splitPart(string s)
{
    pair<string, string> ans = { "", "" };
    size_t i;
    for (i = 0; i < s.size() && s[i] != ' '; ++i) {
        ans.first += s[i];
    }
    for (i += 1; i < s.size(); ++i) {
        ans.second += s[i];
    }
    return ans;
}

double calculate_Compress_rationality(const string& original, const vector<Node>& compressed) {
    size_t original_sz = original.size();
    size_t compressed_sz = compressed.size() * (sizeof(int64_t) * 2 + sizeof(char));
    return (static_cast<double>(compressed_sz) / original_sz);
}

size_t try_to_compressLZ77(istream& input, ostream& output) {
    vector<Node> crypted;
    stringstream buffer;
    buffer << input.rdbuf();
    string s = buffer.str();

    encodeLZ77(s, crypted);
    char flag_of_compression='1';

    if (calculate_Compress_rationality(s, crypted) - 1.0 < std::numeric_limits<double>::epsilon()) {
            output.write(reinterpret_cast<char*>(&flag_of_compression), sizeof(flag_of_compression));
        for (size_t i = 0; i < crypted.size(); ++i) {
            output.write(reinterpret_cast<const char*>(&crypted[i].offset), sizeof(crypted[i].offset));
            output.write(reinterpret_cast<const char*>(&crypted[i].length), sizeof(crypted[i].length));
            output.write(reinterpret_cast<const char*>(&crypted[i].next), sizeof(crypted[i].next));

        }
          return s.size();
    } else {
        cout << "Compression is not effective. Skipping LZ77 compression." << endl;
        output<<'0'<<s;

        return s.size();
    }
}

void decompressLZ77(istream& input, ostream& output) {
    vector<Node> crypted;
    int64_t a, b;
    char c;
    char flag_of_compression;
    input.read(reinterpret_cast<char*>(&flag_of_compression), sizeof(flag_of_compression));
    if('0'==flag_of_compression){output << input.rdbuf();return;}
    while (input.read(reinterpret_cast<char*>(&a), sizeof(a)) && input.read(reinterpret_cast<char*>(&b), sizeof(b)) && input.read(&c, sizeof(c))) {
        crypted.push_back({ a, b, c });
    }
    string s = decodeLZ77(crypted);
    output << s;
}

void arif_compress_with_lz77(std::istream& fin, std::ostream& fout) {

    std::stringstream lz_compressed_data;
    auto start_arif = std::chrono::high_resolution_clock::now();
        fsz = try_to_compressLZ77(fin, lz_compressed_data);
    auto end_arif = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_arif = end_arif - start_arif;
    auto start_lz77 = std::chrono::high_resolution_clock::now();
        arif_compress(lz_compressed_data, fout);
    auto end_lz77 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_lz77 = end_lz77 - start_lz77;
    std::cout<< (duration_arif + duration_lz77).count() << " s" << std::endl;


}

void arif_decompress_with_lz77(std::istream& fin, std::ostream& fout) {

    std::stringstream decompressed_data;

    auto start_arif = std::chrono::high_resolution_clock::now();
    arif_decompress(fin, decompressed_data);
    auto end_arif = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_arif = end_arif - start_arif;

   
    auto start_lz77 = std::chrono::high_resolution_clock::now();
    decompressLZ77(decompressed_data, fout);
    auto end_lz77 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_lz77 = end_lz77 - start_lz77;
    std::cout<< (duration_arif + duration_lz77).count() << " s" << std::endl;
}

bool test_decompression(std::istream& fin){

try{
    std::stringstream decompressed_data;
    std::stringstream fout;
    arif_decompress(fin, decompressed_data);
    decompressLZ77(decompressed_data, fout);

}catch(const std::exception& expecpt){
  std::cerr << "Error in integrity test: " << expecpt.what() << std::endl;
  return false;
}
return true;
}



void process_file(const string& file) {
    std::ifstream inputFile(file, std::ios::binary);
    if (!inputFile.is_open()) {
        std::cerr << "Unable to open file: " << file << std::endl;
        return;
    }
    if((_l||_t||_d)&&file.substr(file.size()-4,file.size())!=".agz"){
       inputFile.close();
        std::cerr << "File is not arhive(.agz): " << file << std::endl;
        return;
    }
    if(_l){
        inputFile.read(reinterpret_cast<char*>(&fsz), sizeof(fsz));
        uint64_t compfsz=filesystem::file_size(file);
        std::cout<<fsz<<" "<<compfsz<<" "<<((double)(fsz-compfsz)*100/fsz)<<"% "<<file.substr(0,file.size()-4)<<std::endl;
    }
   else if (_t) {
        if (test_decompression(inputFile)) {
            std::cout << "Integrity test passed for file: " << file << std::endl;
        } else {
            std::cerr << "Integrity test failed for file: " << file << std::endl;
        }
    } else {
        string outfile;
        if(_d){
            outfile=file.substr(0,file.size()-4);
        } else outfile=file + ".agz";

        std::ofstream outputFile(outfile, std::ios::binary);
        if (!outputFile.is_open()) {
            std::cerr << "Unable to open output file: " << file + ".agz" << std::endl;
            return;
        }

        if (_d) {
            if(!_c)
                arif_decompress_with_lz77(inputFile, outputFile);
            else arif_decompress_with_lz77(inputFile, std::cout);
        } else {
            if(!_c)
                arif_compress_with_lz77(inputFile, outputFile);
            else arif_compress_with_lz77(inputFile, std::cout);
        }


        outputFile.close();
        if(!_k&&!_t&&!_l){
                if (remove(file.c_str())) {
                    std::cerr << "Error deleting file: " << file << std::endl;
                }
            }
    }
    inputFile.close();
}

void process_directory(const std::filesystem::path& directory) {
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_directory()) {
            process_directory(entry.path());
        } else {
            process_file(entry.path().string());
        }
    }
}



int main(int argc, char* argv[]) {
    int opt;
    setlocale(LC_ALL, "Russian");
    while ((opt = getopt(argc, argv, "cdklrt19")) != -1) {
        switch (opt) {
            case '1':
                _1 = true;
                break;
            case '9':
                _9 = true;
                break;
            case 'd':
                _d = true;
                break;
            case 't':
                _t = true;
                break;
            case 'r':
                _r = true;
                break;
            case 'l':
                _l = true;
                break;
            case 'k':
                _k = true;
                break;
            case 'c':
                _c = true;
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " [-c] [-d] [-k] [-l] [-r] [-t] [-1] [-9] [file ...]" << std::endl;
                return 1;
        }
    }

    if (_d && (_1 || _9)) {
        std::cerr << "Error: -d cannot be used with -1 or -9." << std::endl;
        return 1;
    }
    if (_d && _l) {
        std::cerr << "Error: -d cannot be used with -l." << std::endl;
        return 1;
    }
    if (_d && _t) {
        std::cerr << "Error: -d cannot be used with -t." << std::endl;
        return 1;
    }
    if (_l && (_1 || _9)) {
        std::cerr << "Error: -l cannot be used with -1 or -9." << std::endl;
        return 1;
    }
    if (_l && _t) {
        std::cerr << "Error: -l cannot be used with -t." << std::endl;
        return 1;
    }
    if (_t && (_1 || _9)) {
        std::cerr << "Error: -t cannot be used with -1 or -9." << std::endl;
        return 1;
    }
    if(_1 &&_9){
        std::cerr << "Error: -1 cannot be used with -9" << std::endl;
        return 1;
    }

    if(_l){
        cout<<"uncompressed compressed ratio uncompressed_name"<<std::endl;
    }

    std::vector<std::string> files(argv + optind, argv + argc);

    if (files.empty()) {
        std::cerr << "No input files specified." << std::endl;
        return 1;
    }

    for (const auto& file : files) {
    if (_r) {
        struct stat path_stat;
        stat(file.c_str(), &path_stat);
        if (S_ISDIR(path_stat.st_mode)) {
            process_directory(file);
        } else {
            process_file(file);
        }
    } else {
        if (file == "-") {
            if (_t) {
                if (test_decompression(std::cin)) {
                    std::cout << "Integrity test passed." << std::endl;
                } else {
                    std::cerr << "Integrity test failed." << std::endl;
                }
            } else {
                std::ofstream outputFile;
                string outfile="defualt_output.txt";
                if(!_d){outfile+=".agz";}
                if(!_c){
                    outputFile.open(outfile, std::ios::binary);
                }
                if (_d) {
                    if(_c){
                        arif_decompress_with_lz77(std::cin, std::cout);
                    }
                    else  {arif_decompress_with_lz77(std::cin,outputFile);}
                } else {
                    if(_c)
                        arif_compress_with_lz77(std::cin, std::cout);
                    else  arif_compress_with_lz77(std::cin,outputFile);
                }
                if(!_c){
                    outputFile.close();
                }
            }
        } else {
            process_file(file);
        }
    }
}


    return 0;
}
