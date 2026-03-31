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
        :Length(Length), N(N), wx(wx), wy(wy), A(A), t(0.f), height_scale(height_scale),
        h0(N*N, {0.f, 0.f}), ht(N*N, {0.f, 0.f}), Gx(N*N), Gz(N*N), Dx(N*N), Dz(N*N), choppy_lambda(choppy),
        ht_real(N*N), Gx_real(N*N), Gz_real(N*N), Dx_real(N*N), Dz_real(N*N), xyz_coord(N*N*3), gxyz_coord(N*N*3)
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
    void update_realvalue_vector();

    void set_choppy_coefficient(const float lambda) { 
        this->choppy_lambda = lambda;
    }
    const float get_choppy_coefficient() const { 
        return this->choppy_lambda; 
    }
    const std::vector<float> &get_ht_real() const {
        return this->ht_real;
    }
    const std::vector<float> &get_Gx_real() const {
        return this->Gx_real;
    }
    const std::vector<float> &get_Gz_real() const {
        return this->Gz_real;
    }
    const std::vector<float> &get_Dx_real() const {
        return this->Dx_real;
    }
    const std::vector<float> &get_Dz_real() const {
        return this->Dz_real;
    }
    const float get_time() const {
        return this->t;
    }
    uintptr_t get_ht_real_ptr() {  return reinterpret_cast<uintptr_t>(this->ht_real.data());    }
    float* get_Gx_real_ptr() {  return this->Gx_real.data();    }
    float* get_Gz_real_ptr() {  return this->Gz_real.data();    }
    float* get_Dx_real_ptr() {  return this->Dx_real.data();    }
    float* get_Dz_real_ptr() {  return this->Dz_real.data();    }
    uintptr_t get_xyz_ptr() { return reinterpret_cast<uintptr_t>(this->xyz_coord.data()); }
    uintptr_t get_gxyz_ptr() { return reinterpret_cast<uintptr_t>(this->gxyz_coord.data()); }
    inline 
    float get_height_scale(void) const {
        return this->height_scale;
    }
    inline 
    void set_height_scale(const float new_value) {
        this->height_scale = new_value;
    }


private:
    inline 
    void copy_real_value(const ComplexContainer &src_vector, std::vector<float> &dest_vector) {
        for(size_t i = 0; i < src_vector.size(); i++) {
            dest_vector[i] = src_vector[i].real();
        }
    }
    inline 
    void copy_imag_value(const ComplexContainer &src_vector, std::vector<float> &dest_vector) {
        for(size_t i = 0; i < src_vector.size(); i++) {
            dest_vector[i] = src_vector[i].imag();
        }
    }

    inline 
    void form_xyz_array(void) {
        bool calculate_choppy = this->choppy_lambda > 1e-6f ? true : false;
        for(size_t i = 0; i < this->N; i++) {
            for(size_t j = 0; j < this->N; j++) {
                size_t idx = i * N + j;
                float x = ((float)i / (float)(N-1)) * Length - Length / 2.0f;
                float z = ((float)j / (float)(N-1)) * Length - Length / 2.0f;
                float y = this->ht[idx].real();

                this->xyz_coord[idx*3+0] = x;
                this->xyz_coord[idx*3+1] = y;
                this->xyz_coord[idx*3+2] = z;

                if (calculate_choppy) {
                    float choppy_x = this->Dx[idx].real();
                    float choppy_z = this->Dz[idx].real();
                    this->xyz_coord[idx*3+0] += this->choppy_lambda * choppy_x;
                    this->xyz_coord[idx*3+2] += this->choppy_lambda * choppy_z;
                } 


                float dx = this->Gx[idx].real();
                float dz = this->Gz[idx].real();
                float grad_norm = sqrt(dx*dx + 1.f*1.f * dz*dz);
                this->gxyz_coord[idx*3+0] = -dx / grad_norm;
                this->gxyz_coord[idx*3+1] = 1.f / grad_norm;
                this->gxyz_coord[idx*3+2] = -dz / grad_norm;
            }
        }
    }
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
    float t;
    float choppy_lambda;
    float height_scale;
    // Containers to store only the real values, without imagenary part.
    std::vector<float> ht_real;
    std::vector<float> Gx_real;
    std::vector<float> Gz_real;
    std::vector<float> Dx_real;
    std::vector<float> Dz_real;
    std::vector<float> xyz_coord;
    std::vector<float> gxyz_coord;
};

}   // namespace FFTOcean