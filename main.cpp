#include <cstdio>
#include <iostream>
#include <set>
#include <fstream>
#include <map>
#include <ctime>
#include <vector>

#include "htslib/sam.h"
#include "htslib/faidx.h"

using namespace std;
const int WINDOW_SIZE = 5000;
const char *fasta = "result.fa";

struct gc_correction_window_t {
    string chr;
    int pos;
    double gc_percentage;
    mutable int total_read_depth;

    gc_correction_window_t(string chr0, const int pos0)
            :
            chr(std::move(chr0)),
            pos(pos0),
            total_read_depth(0) {
    }

    bool operator<(const gc_correction_window_t &obj) const {
        if (this->chr != obj.chr) return this->chr < obj.chr;
        return this->pos < obj.pos;
    }
};

struct singly_unique_nucleotide {
    string chr;
    char nucleotide;
    int idx;

    singly_unique_nucleotide(string chr0, const int index0, const char nucleotide) :
            chr(std::move(chr0)), idx(index0), nucleotide(nucleotide) {}

    singly_unique_nucleotide(string line) {
        size_t s_pos = 0;
        size_t e_pos = line.find(':');
        chr = line.substr(s_pos, e_pos - s_pos);
        s_pos = e_pos + 1;
        e_pos = line.find(':', e_pos + 1);
        idx = stoi(line.substr(s_pos, e_pos - s_pos));
        nucleotide = line[e_pos + 1];
    }

    bool operator<(const singly_unique_nucleotide &location) const {
        if (this->chr != location.chr) return this->chr < location.chr;
        if (this->idx != location.idx) return this->idx < location.idx;
        return false;
    }

    bool operator==(const singly_unique_nucleotide &location) const {
        return this->chr == location.chr && this->idx == location.idx;
    }

    bool operator>(const singly_unique_nucleotide &location) const {
        return !(*this == location || *this < location);
    }

};

map<singly_unique_nucleotide, int> sun_read(const string &filename);

double get_gc_percent(const char *fasta, const char *region);

int main() {
    map<singly_unique_nucleotide, int> sun_set = sun_read("/Users/emre/CLionProjects/untitled/cmake-build-debug/test");
    //cout << sun_set.size() << endl;
    map<char *, vector<gc_correction_window_t>> read_depth_table;
    set<gc_correction_window_t> active_windows;
    samFile *fp = hts_open("/Volumes/Seagate Expansion Drive/NA12890.recal.bam", "r");
    bam_hdr_t *bam_hdr = sam_hdr_read(fp);
    bam1_t *aln = bam_init1();

    long long unsigned int s = 0;
    clock_t begin = clock();
    string current_chr = "chr1";

    while (sam_read1(fp, bam_hdr, aln) > 0) {
        s++;
        int32_t pos = aln->core.pos + 1; //left most position of alignment in zero based coordianate (+1)
        string chr(bam_hdr->target_name[aln->core.tid]); //contig name (chromosome)
        auto len = static_cast<uint32_t>(aln->core.l_qseq); //length of the read.

        uint8_t *q = bam_get_seq(aln); //quality string
        //uint32_t q2 = aln->core.qual ; //mapping quality

        char *qseq = (char *) malloc(len);
        if (chr.find("chr") == string::npos) {
            chr = "chr" + chr;
        }

        if (current_chr != chr) {
            current_chr = chr;
            for (auto ptr : active_windows) {
                string region = ptr.chr + ":" + to_string(ptr.pos) + "-" + to_string(ptr.pos + WINDOW_SIZE - 1);
                ptr.gc_percentage = get_gc_percent(fasta, region.c_str());
                /*auto key = new char[5];
                sprintf(key, "%.1f", ptr.gc_percentage);
                read_depth_table[key].push_back(ptr);
                free(key);*/
                cout << region << "\t" << ptr.gc_percentage << "\t"
                     << static_cast<double>(ptr.total_read_depth) / WINDOW_SIZE << endl;
            }
            active_windows.clear();
        }

        /*if (flag) {
            cout << ">" << current_chr << " dna:chromosome " << current_chr << ":" << 1 << ":" << region.size() << "\n";
            for (size_t i = 0; i < region.size(); i += 60) {
                cout << region.substr(i, 60) << endl;
            }
            current_chr = chr;
            region.clear();
            read_depth.clear();
        }

        for(int i=0; i<len; i++) {
            qseq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
            while (region.size()<pos+i){
                region.push_back('N');
                read_depth.push_back(0);
            }
            if(region.size()==pos+i){
                region.push_back(qseq[i]);
                read_depth.push_back(1);
            }
            else {
                if (region[(pos+i)] == 'N' || region[(pos+i)] == qseq[i]){
                    region[(pos+i)] = qseq[i];
                    read_depth[(pos+i)]++;
                }
                else if (region[(pos+i)] != qseq[i]) {
                    if(read_depth[(pos+i)] == 0) {
                        region[(pos+i)] = qseq[i];
                        read_depth[(pos+i)]++;
                    } else {
                        read_depth[(pos+i)]--;
                    }
                }
            }
        }*/

        gc_correction_window_t tmp(chr, WINDOW_SIZE * (pos / WINDOW_SIZE));
        auto related_window = active_windows.insert(tmp);

        if (pos + len > related_window.first->pos + WINDOW_SIZE) {
            gc_correction_window_t tmp2(chr, WINDOW_SIZE * (pos / WINDOW_SIZE + 1));
            auto related_window2 = active_windows.insert(tmp2);

            related_window.first->total_read_depth += tmp2.pos - pos;
            related_window2.first->total_read_depth += pos + len - tmp2.pos;
        } else {
            related_window.first->total_read_depth += len;
        }


        singly_unique_nucleotide start_sun(chr, pos, seq_nt16_str[bam_seqi(q, 0)]),
                end_sun(chr, pos + len, seq_nt16_str[bam_seqi(q, len - 1)]);

        if (start_sun > sun_set.rbegin()->first) break;

        auto it = sun_set.lower_bound(start_sun);

        while (it != sun_set.end() && (it->first < end_sun || it->first == end_sun)) {
            if (qseq[it->first.idx - pos] == it->first.nucleotide) {
                it->second++;
            } else {
                // TODO
            }
            it++;
        }

        /*if (s%1000000==0) {
            cout << 100*(double)s/1377710790.18 << endl;
        }*/
        free(qseq);

    }

    clock_t end = clock();
    cout << "Time elapsed " << double(end - begin) / CLOCKS_PER_SEC << endl;
    bam_destroy1(aln);
    hts_close(fp);
    return 0;
}

map<singly_unique_nucleotide, int> sun_read(const string &filename) {
    map<singly_unique_nucleotide, int> result;
    ifstream file(filename);
    if (file.is_open()) {
        string line;
        char delimiter = '\t';
        while (getline(file, line)) {
            size_t s_pos = line.find(delimiter) + 1;
            size_t e_pos = line.find(delimiter, s_pos);
            string token;
            while (e_pos != string::npos) {
                singly_unique_nucleotide sun(line.substr(s_pos, e_pos - s_pos));
                result.insert(make_pair(sun, 0));
                s_pos = e_pos + 1;
                e_pos = line.find(delimiter, e_pos + 1);
            }
        }
    }
    return result;
}

double get_gc_percent(const char *fasta, const char *region) {
    int length;
    double sum = 0;

    faidx_t *fasta_index = fai_load(fasta);
    char *seq = fai_fetch(fasta_index, region, &length);

    if (length > 0) {
        for (int i = 0; i < length; ++i) {
            if (seq[i] == 'G' || seq[i] == 'C') {
                sum++;
            }
        }
    }
    fai_destroy(fasta_index);
    delete[]seq;

    return 100 * (sum / length);
}