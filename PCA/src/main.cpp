#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <cerrno>

#include "cblas.h"

#include "string_utils.h"

// TODO(Colin): Clean up code and use Eigen or something instead of
// raw BLAS calls...

// TODO(Colin): Refactor code so BRDF related stuff goes into one lib,
// and generic PCA stuff is elsewhere...

// Number of BRDF measurements in a single color channel.
static const int NUMEL_1_BRDF_CHANNEL = 90 * 90 * 180;

// Functions operating on arbitrary dimensional vectors ---------------------------------

// Dot product of two N dimensional vectors.
double dot(double* a, double* b, int N){
#ifdef DONT_USE_BLAS
    double accum = 0.0f;
    for (int i = 0; i < N; ++i) {
        accum += (a[i] * b[i]);
    }
    return accum / static_cast<double>(N);
#else
    DDOT(N, a, 1, b, 1);
#endif
}

// Vector addition.  Perform a += b.
void add(double* a, double* b, int N) {
#ifdef DONT_USE_BLAS
    for (int i = 0; i < N; ++i) {
        a[i] += b[i];
    }
#else
    DAXPY(N, 1.0, b, 1, a, 1);
#endif
}

// Vector subtraction.  Perform a -= b.
void sub(double* a, double* b, int N) {
#ifdef DONT_USE_BLAS
    for (int i = 0; i < N; ++i) {
        a[i] -= b[i];
    }
#else
    DAXPY(N, -1.0, b, 1, a, 1);
#endif
}

// Scalar multiplication.  Perform a *= c.
void scalar_mult(double* a, double c, int N) {
#ifndef DONT_USE_BLAS
    for (int i = 0; i < N; ++i) {
        a[i] *= c;
    }
#else
    DSCAL(N, c, 1, a);
#endif
}

// Take the log base e of all elements in a vector.
void component_wise_log(double* a, int N, bool take_log) {
    if (!take_log) {
        return;
    }

    // TODO How to deal with slighly negative values corrupting results?
    const double SMALL_VAL = 1e-10;
    for (int i = 0; i < N; ++i) {

        if (fabs(a[i] < SMALL_VAL)) {
            a[i] = 0.0;
        } else {
            if (a[i] > 0.0) {
                a[i] = log(a[i]);
            } else {
                a[i] = -log(-a[i]);
            }
        }
    }
    assert(errno == 0); // Make sure no cmath errors...
}

// Covariance matrix class ----------------------------------------------------

// Simple covariance matrix class.  Only for small, dense matrices.
class CovMat {
public:
    CovMat(int n_rows, double init_val) : N(n_rows) {
        assert(N > 0);
        data = new double[N * N];
        for (int i = 0; i < N * N; ++i) {
            data[i] = init_val;
        }
    }
    virtual ~CovMat() {
        delete[] data;
        data = NULL;
    }

    double& operator()(int row, int col) {
        assert(row >= 0 && row < N);
        assert(col >= 0 && col < N);
        return data[row * N + col];
    }

    inline int getNumRows() const { return N; }
    inline int getNumCols() const { return N; }
private:
    int N; // matrix dimensions(symmetric matrix).
    double* data;
};

// Misc utilities -------------------------------------------------------------
inline bool IsNaN(double val) {
    return val != val; // Only true for NaN.
}

bool HasANaN(double* buf, int N) {
    for (int i = 0; i < N; ++i) {
        if (IsNaN(buf[i])) {
            return true;
        }
    }
    return false;
}

// Entry point ----------------------------------------------------------------
int main(int argc, char** argv){
    if (argc < 3) {
        std::cout << "This binary computes a covariance matrix from a large set of "
                  << "measured BRDFs. This is the first step in PCA.\n";
        std::cout << "\n";
                  << "Each BRDF is considered as vector from R^(90*180*180), as \n   "
                  << "described in \"A Data Driven Refelctance Model,\" by Matusik "
                  << "Matusik et al.\n"
        std::cout << "\n";
        std::cout << "To actually find the principal components, you still \n"
                  << "need to perform eigenanalysis on this resulting covariance \n"
                  << "matrix.  This is implemented in an included octave/matlab \n"
                  << "script \"pca.m\".\n";
        std::cout << "\n";
        std::cout << "Usage: " << std::endl << argv[0]
                  << " in_file_1, ... , in_file_N, out_name" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "\t--take_natural_log   (boolean) Should we take the log of each BRDF value?\n";
        std::cout << "\t--whiten_data        (boolean) Should we whiten the matrix(subtract mean vector)?\n";
        std::cout << "\t--scale_covariances  (boolean) Should cov(x,y) be dot(x,y) or (1/N)dot(x,y)?\n";
        std::cout << "\t--whiten_before_log  (boolean) Should we whiten raw BRDF or log-BRDF?\n";
        std::cout << "\t--num_color_channels (int)     How many color channels in the BRDF(default to 3)?\n";
        return 1;
    }

    // Parse optional command line flags.
    int num_color_channels = 3;
    // Note if we take the log we MUST do the sqrt later on in the pca.m script!
    bool take_log = true;
    bool whiten_data = true;
    bool scale_covariances = false;
    bool whiten_before_log = true;

    int curr_arg_index = 1;
    for (; curr_arg_index < argc; ++curr_arg_index) {
        std::string curr_arg_str(argv[curr_arg_index]);
        bool parse_ok = true;
        if (StringUtils::startsWith(curr_arg_str, "--take_natural_log")) {
            take_log = StringUtils::parseArg<bool>(curr_arg_str, "--take_natural_log", &parse_ok);
        } else if (StringUtils::startsWith(curr_arg_str, "--whiten_data")) {
            whiten_data = StringUtils::parseArg<bool>(curr_arg_str, "--whiten_data", &parse_ok);
        } else if (StringUtils::startsWith(curr_arg_str, "--scale_covariances")) {
            scale_covariances = StringUtils::parseArg<bool>(curr_arg_str, "--scale_covariances", &parse_ok);
        } else if (StringUtils::startsWith(curr_arg_str, "--whiten_before_log")) {
            whiten_before_log = StringUtils::parseArg<bool>(curr_arg_str, "--whiten_before_log", &parse_ok);
        } else if (StringUtils::startsWith(curr_arg_str, "--num_channels")) {
            num_color_channels = StringUtils::parseArg<int>(curr_arg_str, "--num_channels", &parse_ok);
        } else {
            // There are no optional flags left to read.
            break;
        }
        if (!parse_ok) {
            std::cerr << "Could not parse argument: " << curr_arg_str << std::endl;
            return 1;
        }
    }
    int arg_filenames_start = curr_arg_index;

    // Open output file.
    std::string out_name(argv[argc-1]);
    std::fstream out_file(out_name.c_str(), std::fstream::out);
    std::cout << "Writing to output file: \"" << out_name << "\"." << std::endl;
    const char* WARN_ENDINGS[3] = {".binary", ".brdf", ".sbrdf"};
    bool bad_name = false;
    for (int i = 0; i < 3; ++i) {
        if (StringUtils::endsWith(out_name, WARN_ENDINGS[i])) {
            bad_name = true;
            break;
        }
    }
    if (bad_name) {
        std::cout << "WARNING - The output file name, \"" << out_name << "\" looks like a BRDF file," << std::endl;
        char answer = ' ';
        while (answer != 'y' && answer != 'n') {
            std::cout << "Proceed? (y/n): ";
            std::cin >> answer;
            std::cout << std::endl;
        }
        if (answer == 'n') {
            return 1;
        }
    }

    // Open an series of input files.  Each input file is one row in the data matrix.
    std::cout << "Opening input files:" << std::endl;
    std::vector<std::fstream*> in_files;
    for (int i = arg_filenames_start; i < argc - 1; ++i) {
        std::fstream* fs = new std::fstream(argv[i], std::fstream::binary | std::fstream::in);
        if (!fs->good()) {
            std::cerr << "Could not open file: \"" << argv[i] << "\"" << std::endl;
            fs->close();
            delete fs;
            out_file.close();
            for (size_t j = 0; j < in_files.size(); ++j) {
                in_files[j]->close();
                delete in_files[j];
            }
            return 1;
        }
        in_files.push_back(fs);
        std::cout << "\tOpened input file: \"" << argv[i] << "\"." << std::endl;
    }
    std::cout << "Opened " << in_files.size() << " BRDFs." << std::endl;

    // Declare vars used throughout the rest of the program.
    // Print a bit of status info that the user can look at while things run.
    const int HEADER_SIZE_BYTES = 3 * sizeof(int);
    const int matrix_row_size = NUMEL_1_BRDF_CHANNEL * num_color_channels;
    std::cout << "Each BRDF is being considered as a vector from R^" << matrix_row_size << std::endl;
    std::cout << "Settings:" << std::endl;
    std::cout << "\tNum BRDFs          = " << in_files.size()    << std::endl;
    std::cout << "\tTaking natural log = " << take_log           << std::endl;
    std::cout << "\tWhiten data        = " << whiten_data        << std::endl;
    std::cout << "\tWhiten before log  = " << whiten_before_log  << std::endl;
    std::cout << "\tNum color channels = " << num_color_channels << std::endl;
    std::cout << "\tOutput file        = " << out_name           << std::endl;
    std::cout << "\tScaling covariance = " << scale_covariances  << std::endl;

    double* mean_buf = new double[matrix_row_size]; // Buffer for avg. value.
    double* buf_a = new double[matrix_row_size]; // Buffer for dot product.
    double* buf_b = new double[matrix_row_size]; // Buffer for dot product.

    // Zero out all buffers
    for (int i = 0; i < matrix_row_size; ++i) {
        mean_buf[i] = buf_a[i] = buf_b[i] = 0.0;
    }

    // Find the mean BRDF.
    if (whiten_data) {
        std::cout << "Computing average BRDF for whitening..." << std::endl;
        for (size_t i = 0; i < in_files.size(); ++i) {
            // Read current BRDF info buf_a.
            in_files[i]->seekg(HEADER_SIZE_BYTES);
            in_files[i]->read( (char*)(buf_a), matrix_row_size * sizeof(double));
            if (!whiten_before_log) {
                component_wise_log(buf_a, matrix_row_size, take_log);
            }

            // Add mean_buf += buf_a
            add(mean_buf, buf_a, matrix_row_size);
        }
        scalar_mult(mean_buf, 1.0 / static_cast<double>(in_files.size()), matrix_row_size);
        //if (whiten_before_log) {
        //    component_wise_log(mean_buf, matrix_row_size, take_log);
        //}
        std::cout << "Done computing average BRDF." << std::endl;
    } else {
        std::cout << "Skipped data whitening(aka mean BRDF subtraction)." << std::endl;
    }

    // Make sure we are NaN free
    assert(! HasANaN(mean_buf, matrix_row_size));
    assert(! HasANaN(buf_a, matrix_row_size));
    assert(! HasANaN(buf_b, matrix_row_size));


    // Multiply the matrix A * A^T.

    // We load data successivley info buf_b and buf_a and repeatedly compute
    // dot products.  This computes A * A^T.  This is the (scaled) covariance matrix.
    // The covariance matrix is summetric, so for an N by N covariance matrix, we only need to
    // compute (N^2)/2 dot products.
    int num_rows = static_cast<int>(in_files.size());
    CovMat covariance_matrix(num_rows, -999.0);
    int count = 0;
    const int count_max = (num_rows * (num_rows+1)) / 2;
    std::cout << "Computing covariance matrix entries..." << std::endl;;
    const int REQ_PERC_JUMP_FOR_UPDATE = 5; // Update every X percent
    int prev_perc = -REQ_PERC_JUMP_FOR_UPDATE * 2; // Always update at 0 percent
    for (int r = 0; r < num_rows; ++r) {

        // Read in BRDF buffer for this row.
        in_files[r]->seekg(HEADER_SIZE_BYTES);
        in_files[r]->read( (char*)(buf_a), matrix_row_size * sizeof(double));
        if (whiten_data) {
            if (whiten_before_log) {
                sub(buf_a, mean_buf, matrix_row_size);
                component_wise_log(buf_a, matrix_row_size, take_log);
            } else {
                component_wise_log(buf_a, matrix_row_size, take_log);
                sub(buf_a, mean_buf, matrix_row_size);
            }
        } else {
            component_wise_log(buf_a, matrix_row_size, take_log);
        }

        for (int c = 0; c <= r; ++c) {

            // Read in BRDF buffer for this column.
            in_files[c]->seekg(HEADER_SIZE_BYTES);
            in_files[c]->read( (char*)(buf_b), matrix_row_size * sizeof(double));
            if (whiten_data) {
                if (whiten_before_log) {
                    sub(buf_b, mean_buf, matrix_row_size);
                    component_wise_log(buf_b, matrix_row_size, take_log);
                } else {
                    component_wise_log(buf_b, matrix_row_size, take_log);
                    sub(buf_b, mean_buf, matrix_row_size);
                }
            } else {
                component_wise_log(buf_b, matrix_row_size, take_log);
            }

            // Compute the large dot product.
            double value = dot(buf_a, buf_b, matrix_row_size) ;
            if (scale_covariances) {
                value /= static_cast<double>(matrix_row_size);
            }
            assert(!IsNaN(value));

            covariance_matrix(r, c) = value;
            covariance_matrix(c, r) = value;

            // Write to the output file, and log our progress to stdout.
            int perc = static_cast<int>((static_cast<double>(count) / static_cast<double>(count_max)) * 100.0);
            if (perc >= prev_perc + REQ_PERC_JUMP_FOR_UPDATE) {
                std::cout << "\tThe computation is: " << perc << " percent complete." << std::endl;
                std::cout << "\t\tMost recent covariance entry: cov(" << r << ", " << c << ") = " << value << std::endl;
                prev_perc = perc;
            }
            ++count;
        }
    }
    std::cout << "Done computing covariance matrix. Outputting results to file: " << out_name << std::endl;
    // Write results to file.
    for (int r = 0; r < covariance_matrix.getNumRows(); ++r) {
        for (int c = 0; c < covariance_matrix.getNumCols(); ++c) {
            out_file << covariance_matrix(r,c) << "    ";
        }
        out_file << std::endl;
    }

    std::cout << "All done.  Results were written to: " << out_name << std::endl;

    // Cleanup mem.
    delete[] buf_a;
    delete[] buf_b;
    delete[] mean_buf;

    // Close all the streams.
    out_file.close();
    for (size_t i = 0; i < in_files.size(); ++i) {
        in_files[i]->close();
        delete in_files[i];
    }
    out_file.close();

    return 0;
}
