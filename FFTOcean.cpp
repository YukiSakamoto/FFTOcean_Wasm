#include "pocketfft_hdronly.h"
#include <vector>
#include <complex>
#include <random>
#include "FFTOcean.hpp"

namespace FFTOceanNamespace {

std::default_random_engine generator;
std::normal_distribution<float> dist(0.0, 1.0);

void FFTOcean::calculate_h0(void)
{
    for(int i = 0; i < this->N; i++) {
        for(int j = 0; j < this->N; j++) {
            float kx = 2.f * M_PI * (i-N/2.f) / Length;
            float ky = 2.f * M_PI * (j-N/2.f) / Length;
            float p = this->phillips(kx, ky);
            float r1 = dist(generator);
            float r2 = dist(generator);
            std::complex<float> c = sqrtf(p*0.5f) * std::complex<float>(r1,r2);
            this->h0[i*N + j] = c;
        }
    }
}

void FFTOcean::calculate_ht(float t) 
{
    for(int i = 0; i < this->N; i++) {
        for(int j = 0; j < this->N; j++) {
            int idx = i*this->N + j;
            int i_inv = (N-i) % N;
            int j_inv = (N-j) % N;
            int idx_inv = i_inv * N + j_inv;

            float kx = 2.f * M_PI * (i-N/2.f) / Length;
            float ky = 2.f * M_PI * (j-N/2.f) / Length;
            float k = sqrtf(kx*kx + ky*ky);

            float omega = sqrtf(G*k);
            std::complex<float> phase(cosf(omega*t), sinf(omega*t));
            ht[idx] = h0[idx] * phase + std::conj(h0[idx_inv]) * std::conj(phase);
        }
    }
}

void FFTOcean::calculate_grad(void)
{
    bool calculate_choppy = this->choppy_lambda > 1e-6f ? true : false;

    const int n = this->N;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            int idx = i * n + j;
            float kx = 2.f * M_PI * (i - n/2.f) / this->Length;
            float kz = 2.f * M_PI * (j - n/2.f) / this->Length;
            float k = sqrtf(kx*kx + kz*kz);

            std::complex<float> i_k(0, 1);
            this->Gx[idx] = ht[idx] * (i_k * kx);
            this->Gz[idx] = ht[idx] * (i_k * kz);

            if (calculate_choppy) {
                if (1e-6f < k) {
                    this->Dx[idx] = ht[idx] * (i_k * (kx / k));
                    this->Dz[idx] = ht[idx] * (i_k * (kz / k));
                } else {
                    this->Dx[idx] = this->Dz[idx] = 0.f;
                }
            }
        }
    }
}


void FFTOcean::iFFT() {
    pocketfft::shape_t shape = {N, N};
    pocketfft::stride_t stride = {(long)(N*sizeof(std::complex<float>)), sizeof(std::complex<float>)};
    pocketfft::shape_t axes = {0, 1};
    const float nn_inv = 1.f/(N*N);
    pocketfft::c2c(shape, stride, stride, axes, false, ht.data(), ht.data(), nn_inv * this->height_scale);
    pocketfft::c2c(shape, stride, stride, axes, false, Gx.data(), Gx.data(), nn_inv * this->height_scale);
    pocketfft::c2c(shape, stride, stride, axes, false, Gz.data(), Gz.data(), nn_inv * this->height_scale);

    bool calculate_choppy = this->choppy_lambda > 1e-6f ? true : false;
    if (calculate_choppy) {
        pocketfft::c2c(shape, stride, stride, axes, false, Dx.data(), Dx.data(), nn_inv);
        pocketfft::c2c(shape, stride, stride, axes, false, Dz.data(), Dz.data(), nn_inv);
    }
}

void FFTOcean::processSign() {
    for(int i = 0; i < this->N; i++) {
        for(int j = 0; j < this->N; j++) {
            size_t idx = i*N + j;
            float sign = ((i+j) % 2 == 0) ? 1.f : -1.f;
            this->ht[idx] *= sign;
            this->Gx[idx] *= sign;
            this->Gz[idx] *= sign;
            this->Dx[idx] *= sign;
            this->Dz[idx] *= sign;
        }
    }
}

void FFTOcean::update_realvalue_vector()
{
    bool calculate_choppy = this->choppy_lambda > 1e-6f ? true : false;
    this->copy_real_value(ht, ht_real);
    this->copy_real_value(Gx, Gx_real);
    this->copy_real_value(Gz, Gz_real);
    if (calculate_choppy) {
        this->copy_real_value(Dx, Dx_real);
        this->copy_real_value(Dz, Dz_real);
    }
}

void FFTOcean::Update(float t)
{
    this->t = t;
    this->calculate_ht(t);
    this->calculate_grad();
    this->iFFT();
    this->processSign();
    this->update_realvalue_vector();
    this->form_xyz_array();
}

}   // namespace FFTOcean