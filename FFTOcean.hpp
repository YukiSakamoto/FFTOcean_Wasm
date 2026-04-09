#include "pocketfft_hdronly.h"
#include <vector>
#include <complex>
#include <random>

namespace FFTOceanNamespace
{

const float G = 9.81f;  // gravity constant

class FFTOcean
{
public:
    typedef std::vector<std::complex<float>> ComplexContainer;
    FFTOcean(size_t Length, size_t N, float wx = 0.f, float wy = 0.f, float A = 0.0005f, float choppy = 0.f, float height_scale = 1.0)
        :Length(Length), N(N), wx(wx), wy(wy), A(A), t(0.f), height_scale(height_scale), choppy_lambda(choppy),
        h0(N*N, {0.f, 0.f}), ht(N*N, {0.f, 0.f}), Gx(N*N), Gz(N*N), 
        Dx(N*N), Dz(N*N), Dxx(N*N), Dxz(N*N), Dzx(N*N), Dzz(N*N),
        xyz_coord(N*N*3), gxyz_coord(N*N*3)
    {
        calculate_h0();
        Update(0.f);
    }
    inline
    float phillips(const float kx, const float ky) {
        float k_sq = kx*kx + ky*ky;
        if (k_sq < 1.e-8) { return 0.f; }
        float k_norm = sqrt(k_sq);
        float w_norm = sqrt(wx*wx + wy*wy);
        if (w_norm < 1.e-6f) {
            return 0.f;
        }
        float L_max = (w_norm*w_norm) / G;
        float cos_factor = (kx/k_norm) * (wx/w_norm) + (ky/k_norm) * (wy/w_norm);
        float ret = A * exp(-1.f / (k_sq*L_max*L_max)) / (k_sq*k_sq) * cos_factor*cos_factor;
        return ret;
    }
    inline
    float operator()(const size_t i, const size_t j) const {
        return ht[i*N + j].real();
    }
    void calculate_h0(void);
    void calculate_ht(float t);
    void calculate_grad();
    void Update(float t);
    void iFFT();
    void processSign();

    void set_choppy_coefficient(const float lambda) { 
        this->choppy_lambda = lambda;
    }
    const float get_choppy_coefficient() const { 
        return this->choppy_lambda; 
    }
    const float get_time() const {
        return this->t;
    }

    inline 
    float get_height_scale(void) const {
        return this->height_scale;
    }
    inline 
    void set_height_scale(const float new_value) {
        this->height_scale = new_value;
    }

    uintptr_t get_xyz_ptr() { return reinterpret_cast<uintptr_t>(this->xyz_coord.data()); }
    uintptr_t get_gxyz_ptr() { return reinterpret_cast<uintptr_t>(this->gxyz_coord.data()); }
    void form_xyz_array();

private:
    size_t Length;
    size_t N;
    float A;
    float wx, wy;
    // Containers fr wave heights
    ComplexContainer h0;
    ComplexContainer ht;
    // Containers for Gradients
    ComplexContainer Gx;
    ComplexContainer Gz;
    // Container for choppy waves
    ComplexContainer Dx;
    ComplexContainer Dz;
    ComplexContainer Dxx, Dxz, Dzx, Dzz;
    float t;
    float choppy_lambda;
    float height_scale;
    // Containers to store only the real values, without imagenary part.
    std::vector<float> xyz_coord;
    std::vector<float> gxyz_coord;
};

}   // namespace FFTOcean