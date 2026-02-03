#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <memory>
#include <algorithm>
#include <ctime>
#include <map>
#include <sstream>
#include <zlib.h>
#include <chrono>

#if defined(_MSC_VER) && _MSC_VER < 1914
    #include <experimental/filesystem>
    namespace fs = std::experimental::filesystem;
#else
    #include <filesystem>
    namespace fs = std::filesystem;
#endif

// Constants for the sorter
constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer for reading
constexpr size_t MAX_MEMORY = 4ULL * 1024 * 1024 * 1024;  // 4GB chunk size

struct ProgramParams {
    std::vector<std::string> reads_files;  // Multiple forward reads
    std::vector<std::string> pairs_files;  // Multiple reverse reads (optional)
    std::string output1_file;    // Output for forward reads
    std::string output2_file;    // Output for reverse reads
    size_t barcode_start;        // Position where barcode starts (0-based)
    size_t barcode_length;       // Length of the barcode
    std::string temp_dir;        // Directory for temporary files
    
    ProgramParams() : 
        output1_file("output1.fastq"),
        output2_file("output2.fastq"),
        barcode_start(3),
        barcode_length(15),
        temp_dir(fs::temp_directory_path().string()) {}
};

struct FastqEntry {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
    
    std::string getBarcode(size_t start, size_t length) const {
        if (start >= sequence.size()) return std::string();
        size_t avail = sequence.size() - start;
        return sequence.substr(start, std::min(length, avail));
    }
    
    // Estimate memory usage of this entry
    size_t memoryUsage() const {
        return header.capacity() + sequence.capacity() + 
               plus.capacity() + quality.capacity();
    }
};

// Helper function to check if file is gzipped
bool isGzipped(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) return false;
    
    unsigned char magic[2];
    file.read(reinterpret_cast<char*>(magic), 2);
    return (file && magic[0] == 0x1f && magic[1] == 0x8b);
}

// Helper class for buffered FASTQ reading
class FastqReader {
    gzFile gz_file;
    std::ifstream text_file;
    bool is_gzipped;
    char buffer[BUFFER_SIZE];

public:
    FastqReader(const std::string& filename) {
        is_gzipped = isGzipped(filename);
        if (is_gzipped) {
            gz_file = gzopen(filename.c_str(), "r");
            if (!gz_file) {
                throw std::runtime_error("Cannot open gzipped file: " + filename);
            }
        } else {
            text_file.open(filename);
            if (!text_file) {
                throw std::runtime_error("Cannot open file: " + filename);
            }
        }
    }

    ~FastqReader() {
        if (is_gzipped) {
            gzclose(gz_file);
        }
    }

    bool getline(std::string& line) {
        line.clear();
        if (is_gzipped) {
            char c;
            while (gzread(gz_file, &c, 1) == 1) {
                if (c == '\n') return true;
                line += c;
            }
            return !line.empty();
        } else {
            return static_cast<bool>(std::getline(text_file, line));
        }
    }

    bool getEntry(FastqEntry& entry) {
        if (!getline(entry.header)) return false;
        if (!getline(entry.sequence)) return false;
        if (!getline(entry.plus)) return false;
        if (!getline(entry.quality)) return false;
        return true;
    }
};

// Helper class for FASTQ writing
class FastqWriter {
    std::ofstream file;

public:
    FastqWriter(const std::string& filename) : file(filename) {
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }
    }

    void writeEntry(const FastqEntry& entry) {
        file << entry.header << '\n'
             << entry.sequence << '\n'
             << entry.plus << '\n'
             << entry.quality << '\n';
        if (!file) {
            throw std::runtime_error("Failed to write to FASTQ file");
        }
    }
};

// Class to manage a chunk file
class ChunkFile {
    std::ifstream file;
    FastqEntry current_entry;
    bool has_next = false;
    std::string filename;

public:
    ChunkFile(const std::string& fname) : file(fname), filename(fname) {
        has_next = readNext();
    }

    bool readNext() {
        std::string line;
        if (!std::getline(file, current_entry.header)) {
            has_next = false;
            return false;
        }
        if (!std::getline(file, current_entry.sequence) ||
            !std::getline(file, current_entry.plus) ||
            !std::getline(file, current_entry.quality)) {
            throw std::runtime_error("Incomplete FASTQ entry in chunk file: " + filename);
        }
        has_next = true;
        return true;
    }

    const FastqEntry& current() const { return current_entry; }
    bool hasNext() const { return has_next; }
    const std::string& getFilename() const { return filename; }
};

// Custom comparator for chunk files in priority queue
class ChunkComparator {
    const ProgramParams& params;

public:
    explicit ChunkComparator(const ProgramParams& p) : params(p) {}

    bool operator()(const std::unique_ptr<ChunkFile>& a, 
                   const std::unique_ptr<ChunkFile>& b) {
        return a->current().getBarcode(params.barcode_start, params.barcode_length) > 
               b->current().getBarcode(params.barcode_start, params.barcode_length);
    }
};

struct FastqPair {
    FastqEntry forward;
    FastqEntry reverse;
    
    std::string getBarcode(size_t start, size_t length) const {
        return forward.getBarcode(start, length);
    }
    
    size_t memoryUsage() const {
        return forward.memoryUsage() + reverse.memoryUsage();
    }
};

class FastqSorter {
    std::string temp_dir;
    std::vector<std::string> chunk_files_forward;
    std::vector<std::string> chunk_files_reverse;
    ProgramParams params;
    bool has_pairs;
    size_t total_size;
    size_t processed_size;
    int last_progress;
    std::chrono::steady_clock::time_point start_time;
    std::vector<size_t> reads_per_file;
    size_t total_reads;

    std::string formatDuration(std::chrono::steady_clock::time_point now) {
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
        auto hours = elapsed / 3600;
        auto minutes = (elapsed % 3600) / 60;
        auto seconds = elapsed % 60;
        std::stringstream ss;
        if (hours > 0) ss << hours << "h ";
        if (minutes > 0) ss << minutes << "m ";
        ss << seconds << "s";
        return ss.str();
    }

    void reportProgress(int current_progress) {
        if (current_progress > last_progress) {
            auto now = std::chrono::steady_clock::now();
            std::cout << "Progress: " << current_progress << "% (Time elapsed: " 
                     << formatDuration(now) << ")" << std::endl;
            last_progress = current_progress;
        }
    }

public:
    FastqSorter(const ProgramParams& p) 
        : params(p), has_pairs(!p.pairs_files.empty()), 
          total_size(0), processed_size(0), last_progress(-10),
          start_time(std::chrono::steady_clock::now()) {
        auto timestamp = std::chrono::system_clock::now().time_since_epoch().count();
        temp_dir = params.temp_dir + "/fastq_sort_" + std::to_string(timestamp);
        if (!fs::create_directory(temp_dir)) {
            throw std::runtime_error("Failed to create temporary directory: " + temp_dir);
        }
        
        // Calculate total file size
        for (const auto& file : params.reads_files) {
            total_size += fs::file_size(file);
        }
        for (const auto& file : params.pairs_files) {
            total_size += fs::file_size(file);
        }
    }

    ~FastqSorter() {
        // Clean up temporary files
        for (const auto& file : chunk_files_forward) {
            fs::remove(file);
        }
        for (const auto& file : chunk_files_reverse) {
            fs::remove(file);
        }
        fs::remove(temp_dir);
    }

    // Split input file into sorted chunks
    void createSortedChunks() {
        std::vector<std::unique_ptr<FastqReader>> readers_forward;
        std::vector<std::unique_ptr<FastqReader>> readers_reverse;
        
        // Initialize read counters
        reads_per_file.resize(params.reads_files.size(), 0);
        total_reads = 0;
        
        // Open all forward read files
        for (const auto& file : params.reads_files) {
            readers_forward.push_back(std::make_unique<FastqReader>(file));
        }
        
        // Open all reverse read files if we have pairs
        if (has_pairs) {
            for (const auto& file : params.pairs_files) {
                readers_reverse.push_back(std::make_unique<FastqReader>(file));
            }
        }

        std::vector<FastqPair> chunk;
        size_t current_memory = 0;
        size_t chunk_count = 0;
        
        bool more_reads = true;
        while (more_reads) {
            more_reads = false;
            
            // Try to read from each file
            for (size_t i = 0; i < readers_forward.size(); ++i) {
                FastqPair pair;
                if (readers_forward[i]->getEntry(pair.forward)) {
                    more_reads = true;
                    reads_per_file[i]++;  // Count reads from this file
                    total_reads++;        // Increment total read count
                    
                    if (has_pairs) {
                        if (!readers_reverse[i]->getEntry(pair.reverse)) {
                            throw std::runtime_error("Paired file " + params.pairs_files[i] + 
                                                  " has fewer reads than " + params.reads_files[i]);
                        }
                    }
                    
                    // Update progress
                    processed_size += pair.memoryUsage();
                    int current_progress = static_cast<int>((processed_size * 10) / total_size);
                    // reportProgress(current_progress);
                    
                    chunk.push_back(std::move(pair));
                    current_memory += chunk.back().memoryUsage();
                    
                    if (current_memory >= MAX_MEMORY) {
                        writeSortedChunk(chunk, chunk_count++);
                        chunk.clear();
                        current_memory = 0;
                    }
                }
            }
        }

        // Print read count information
        std::cout << "\nRead counts per file:" << std::endl;
        for (size_t i = 0; i < params.reads_files.size(); ++i) {
            std::cout << params.reads_files[i] << ": " 
                     << reads_per_file[i] << " reads" << std::endl;
        }
        std::cout << "\nTotal reads processed: " << total_reads << std::endl << std::endl;

        if (!chunk.empty()) {
            writeSortedChunk(chunk, chunk_count++);
        }
    }

    void writeSortedChunk(std::vector<FastqPair>& chunk, size_t chunk_num) {
        // Sort chunk by barcode (using forward reads)
        std::sort(chunk.begin(), chunk.end(),
                 [this](const FastqPair& a, const FastqPair& b) {
                     return a.getBarcode(params.barcode_start, params.barcode_length) < 
                            b.getBarcode(params.barcode_start, params.barcode_length);
                 });

        std::string chunk_file_forward = temp_dir + "/chunk_" + 
                                       std::to_string(chunk_num) + "_1.fastq";
        chunk_files_forward.push_back(chunk_file_forward);
        FastqWriter writer_forward(chunk_file_forward);

        std::unique_ptr<FastqWriter> writer_reverse;
        if (has_pairs) {
            std::string chunk_file_reverse = temp_dir + "/chunk_" + 
                                           std::to_string(chunk_num) + "_2.fastq";
            chunk_files_reverse.push_back(chunk_file_reverse);
            writer_reverse = std::make_unique<FastqWriter>(chunk_file_reverse);
        }

        for (const auto& pair : chunk) {
            writer_forward.writeEntry(pair.forward);
            if (has_pairs) {
                writer_reverse->writeEntry(pair.reverse);
            }
        }
    }

    // Merge sorted chunks into output file
    void mergeSortedChunks() {
        
    FastqWriter writer_forward(params.output1_file);
    std::unique_ptr<FastqWriter> writer_reverse;
    if (has_pairs) {
        writer_reverse = std::make_unique<FastqWriter>(params.output2_file);
    }

    using ChunkFilePtr = std::unique_ptr<ChunkFile>;
    std::priority_queue<ChunkFilePtr,
                        std::vector<ChunkFilePtr>,
                        ChunkComparator> pq((ChunkComparator(params)));

    std::map<ChunkFile*, std::unique_ptr<ChunkFile>> reverse_chunks;
    size_t loaded_chunks = 0;
    const size_t BATCH_SIZE = 240; // should be safe on Windows / Linux / MacOS
    size_t total_chunks = chunk_files_forward.size();
    size_t current_batch_start = 0;
    size_t merged_reads = 0;

    // Open all chunk files
    while (current_batch_start < total_chunks) {
        size_t batch_end = std::min(current_batch_start + BATCH_SIZE, total_chunks);

        // Load a batch of chunks
        for (size_t i = current_batch_start; i < batch_end; ++i) {
            try {
                auto forward_chunk = std::make_unique<ChunkFile>(chunk_files_forward[i]);
                if (has_pairs) {
                    auto reverse_chunk = std::make_unique<ChunkFile>(chunk_files_reverse[i]);
                    if (forward_chunk->hasNext() && reverse_chunk->hasNext()) {
                        reverse_chunks[forward_chunk.get()] = std::move(reverse_chunk);
                        pq.push(std::move(forward_chunk));
                        loaded_chunks++;
                    }
                } else if (forward_chunk->hasNext()) {
                    pq.push(std::move(forward_chunk));
                    loaded_chunks++;
                }
            } catch (const std::exception& e) {
                std::cerr << "Warning: Failed to load chunk " << i << ": " << e.what() << std::endl;
            }
        }

        std::cout << "Successfully loaded " << loaded_chunks << " out of " 
                  << chunk_files_forward.size() << " chunks." << std::endl;

        size_t merged_reads = 0;
        size_t processed_chunks = 0;

        while (!pq.empty()) {
            ChunkFilePtr chunk = std::move(const_cast<ChunkFilePtr&>(pq.top()));
            pq.pop();

            FastqEntry current_forward = chunk->current();

            if (has_pairs) {
                auto& reverse_chunk = reverse_chunks[chunk.get()];
                FastqEntry current_reverse = reverse_chunk->current();

                writer_forward.writeEntry(current_forward);
                writer_reverse->writeEntry(current_reverse);
                merged_reads++;

                bool forward_has_next = chunk->readNext();
                bool reverse_has_next = reverse_chunk->readNext();

                if (forward_has_next != reverse_has_next) {
                    throw std::runtime_error("Chunks out of sync: " + chunk->getFilename() + 
                                             " and its reverse pair");
                }

                if (forward_has_next && reverse_has_next) {
                    pq.push(std::move(chunk));
                } else {
                    processed_chunks++;
                }
            } else {
                writer_forward.writeEntry(current_forward);
                merged_reads++;

                if (chunk->readNext()) {
                    pq.push(std::move(chunk));
                } else {
                    processed_chunks++;
                }
            }

            if (merged_reads % 1000000 == 0) {
                std::cout << "Merged " << merged_reads << " reads so far..." << std::endl;
            }

            int current_progress = 50 + static_cast<int>((processed_chunks * 50) / total_chunks);
            // reportProgress(current_progress);
        }

        reverse_chunks.clear();
        while (!pq.empty()) pq.pop();

        current_batch_start = batch_end;
    }

    // std::cout << "\nTotal reads merged: " << merged_reads << std::endl;
}

    void sort() {
        std::cout << "Creating sorted chunks..." << std::endl;
        createSortedChunks();
        
        std::cout << "Created " << chunk_files_forward.size() << " chunks." << std::endl;
        
        std::cout << "Merging sorted chunks..." << std::endl;
        mergeSortedChunks();
        
        auto end_time = std::chrono::steady_clock::now();
        std::cout << "Sort complete! Total time: " << formatDuration(end_time) << std::endl;
    }
};

ProgramParams parseCommandLine(int argc, char* argv[]) {
    ProgramParams params;
    std::map<std::string, std::string> args;
    
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 >= argc) {
            throw std::runtime_error("Missing value for argument " + std::string(argv[i]));
        }
        args[argv[i]] = argv[i + 1];
    }
    
    // Required parameter
    if (args.find("--in") == args.end()) {
        throw std::runtime_error("Read 1 file is required (--in)");
    }
    
    // Split comma-separated files
    std::stringstream ss(args["--in"]);
    std::string file;
    while (std::getline(ss, file, ',')) {
        if (!file.empty()) {
            params.reads_files.push_back(file);
        }
    }
    
    // Optional paired files
    if (args.find("--in-pairs") != args.end()) {
        std::stringstream ss_pairs(args["--in-pairs"]);
        while (std::getline(ss_pairs, file, ',')) {
            if (!file.empty()) {
                params.pairs_files.push_back(file);
            }
        }
        
        if (params.pairs_files.size() != params.reads_files.size()) {
            throw std::runtime_error("Number of paired files must match number of read files");
        }
    }
    
    // Calculate total size
    size_t total_size = 0;
    for (const auto& file : params.reads_files) {
        total_size += fs::file_size(file);
    }
    for (const auto& file : params.pairs_files) {
        total_size += fs::file_size(file);
    }
    
    if (args.find("--bc-start") != args.end()) {
        params.barcode_start = std::stoul(args["--bc-start"]);
    }
    
    if (args.find("--bc-len") != args.end()) {
        params.barcode_length = std::stoul(args["--bc-len"]);
    }
    
    if (args.find("--out") != args.end()) {
        params.output1_file = args["--out"];
    }
    
    if (args.find("--out-pairs") != args.end()) {
        params.output2_file = args["--out-pairs"];
    }
    
    if (args.find("--temp") != args.end()) {
        params.temp_dir = args["--temp"];
        if (!fs::exists(params.temp_dir)) {
            throw std::runtime_error("Temporary directory does not exist: " + params.temp_dir);
        }
        if (!fs::is_directory(params.temp_dir)) {
            throw std::runtime_error("Specified temp path is not a directory: " + params.temp_dir);
        }
        if (fs::space(params.temp_dir).available < MAX_MEMORY * 2) {  // Ensure minimum space
            throw std::runtime_error("Insufficient space in temporary directory");
        }
    }
    
    return params;
}

int main(int argc, char* argv[]) {
    try {
        auto params = parseCommandLine(argc, argv);
        FastqSorter sorter(params);
        sorter.sort();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Usage: " << argv[0] << " [options]\nOptions:\n"
                 << "  --in str         (e.g., 'input1.fastq' or 'input1.fastq,input2.fastq')\n"
                 << "  --in-pairs str   (e.g., 'input1_rev.fastq' or 'input1_rev.fastq,input2_rev.fastq')\n"
                 << "  --bc-start int   Barcode start position (Zero-indexed, default: 0)\n"
                 << "  --bc-len int     Barcode length (default: 18)\n"
                 << "  --out str        Output file for sorted read 1\n"
                 << "  --out-pairs str  Output file for sorted read 2 (if using --in2)\n"
                 << "  --temp str       Path to temporary directory for storing chunk files (recommend .)\n" << std::endl;
        return 1;
    }

    return 0;
}
