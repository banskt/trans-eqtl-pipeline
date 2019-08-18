#include <algorithm>
#include <cmath>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>


#include "gzstream.h"
#include "readerwriterqueue.h"
#include "atomicops.h"

namespace mc = moodycamel;

template<typename Out> void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim);
std::string header_to_str(std::vector<std::string>&);
bool startsWith(const std::string& haystack, const std::string& needle);
double get_dose(std::stringstream& , size_t);
double PLtoPP(std::string x);
void print_help();

// VCF constants
const size_t REF = 3;
const size_t ALT = 4;
const size_t FORMAT = 8;
const size_t ENTRY_START = FORMAT + 1;
const char ENTRY_SEP = ':';
const char PL_SEP = ',';
const std::string PL = "PP";


int main (int argc, char *argv[]) {


    // taskqueue constants
    const size_t Q_SIZE = 1000; 
    const std::string QUIT_STR = "~";

    if (argc < 3) {
        print_help();
        return 1;
    }

    // initialize input/ouput streams
    igzstream infile(argv[1]);
    ogzstream outfile(argv[2]);

    // create task queue
    mc::BlockingReaderWriterQueue<std::string> lines(Q_SIZE);

    // set up consumer thread
    std::thread reader([&]() {
        std::string line;
        std::string entry;
        std::string fmt_entry;
        std::stringstream ss_line;

        bool decl = false;
        std::vector<std::string> header;

        ss_line.setf(std::ios::fixed,std::ios::floatfield);
        ss_line.precision(2);

        while (true) {
            lines.wait_dequeue(line);

            // output any remaining strings and then quit
            if (line.compare(QUIT_STR) == 0) {
                outfile << ss_line.rdbuf();
                break;
            }

            // buffer the header data
            if (startsWith(line, "##")) {
                if (!decl) {
                    ss_line << line << std::endl;
                    decl = true;
                } else {
                    header.push_back(line);
                }
                continue;
            } else if (startsWith(line, "#")) {
                header.push_back("##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected Posterior Genotype\">");
                header.push_back(line);
                std::sort(header.begin(), header.end());
                ss_line << header_to_str(header);
                continue;
            }

            // grab information for variant
            std::vector<std::string> row = split(line, '\t');

            // buffer the variant meta data
            for (size_t i = 0; i < FORMAT; i++) {
                ss_line << row[i] << '\t';
            }

            // find the PL index in the format string
            size_t pl_idx = 0;
            size_t idx = 0;
            std::stringstream fmt_ss(row[FORMAT]);
            while (std::getline(fmt_ss, fmt_entry, ENTRY_SEP)) {
                if (fmt_entry.compare(PL) == 0) {
                    pl_idx = idx;
                    break;
                }
                idx++;
            }

            // append dosage info to format string
            ss_line << row[FORMAT] << ":DS" << '\t';

            // process likelihoods for each sample and convert to dosage
            for (size_t i = ENTRY_START; i < row.size() - 1; i++) {
                std::stringstream smpl_ss(row[i]);

                double dose = get_dose(smpl_ss, pl_idx);
                ss_line << row[i] << ENTRY_SEP << dose << '\t';
            }
            // last entry needs newline isntead of tab
            size_t i = row.size() - 1;
            std::stringstream smpl_ss(row[i]);
            double dose = get_dose(smpl_ss, pl_idx);

            ss_line << row[i] << ENTRY_SEP << dose << std::endl;
            outfile << ss_line.rdbuf();
            ss_line.str(std::string());
        }
    });

    // read in a bunch of data to process
    // 'sleep' once queue is full
    // TODO: replace with proper condition variable instead of sleep/spool
    std::string line;
    while (std::getline(infile, line)) {
        while(!lines.try_enqueue(line)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }

    // no lines left push quit signal
    while(!lines.try_enqueue(QUIT_STR)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    // join thread and close files
    reader.join();
    outfile.close();
    infile.close();

    return 0;
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

std::string header_to_str(std::vector<std::string>& header) {
    std::stringstream ss_line;
    for (std::string & str : header) {
        ss_line << str << std::endl;
    }
    return ss_line.str();
}

bool startsWith(const std::string& haystack, const std::string& needle) {
        return needle.length() <= haystack.length() 
                    && std::equal(needle.begin(), needle.end(), haystack.begin());
}


double get_dose(std::stringstream& smpl_ss, size_t pl_idx) {
    size_t k = 0;
    double pls [3];
    std::string entry;
    std::string token;

    while (std::getline(smpl_ss, token, ENTRY_SEP)) {
        if (k == pl_idx) {
            std::stringstream ss;
            ss.str(token);
            for (size_t j = 0; j < 3; j++) {
                std::getline(ss, entry, PL_SEP);
                pls[j] = PLtoPP(entry);
            }
            break;
        }
        k++;
    }

    return (pls[1] + 2 * pls[2]) / (pls[0] + pls[1] + pls[2]);
}

double PLtoPP(std::string x) {
    std::string::size_type sz; 
    return std::pow(10.0, std::stod(x, &sz) / -10.0);
}

void print_help() {
    std::cerr << "Usage: add_dosage [OPTIONS] INPUT.VCF.GZ OUTPUT.VCF.GZ" << std::endl;
    std::cerr << "Try 'add_dosage --help' for more options" << std::endl;
    return;
}
