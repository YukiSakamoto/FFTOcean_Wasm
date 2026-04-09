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
                    this->Dxx[idx] = i_k * kx * this->Dx[idx];
                    this->Dxz[idx] = i_k * kz * this->Dx[idx];
                    this->Dzx[idx] = i_k * kx * this->Dz[idx];
                    this->Dzz[idx] = i_k * kz * this->Dz[idx];
                } else {
                    this->Dx[idx] = this->Dz[idx] = 0.f;
                    this->Dxx[idx] = this->Dxz[idx] = this->Dzx[idx] = this->Dzz[idx] = 0.f;
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
        pocketfft::c2c(shape, stride, stride, axes, false, Dxx.data(), Dxx.data(), nn_inv);
        pocketfft::c2c(shape, stride, stride, axes, false, Dxz.data(), Dxz.data(), nn_inv);
        pocketfft::c2c(shape, stride, stride, axes, false, Dzx.data(), Dzx.data(), nn_inv);
        pocketfft::c2c(shape, stride, stride, axes, false, Dzz.data(), Dzz.data(), nn_inv);
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

            this->Dxx[idx] *= sign;
            this->Dxz[idx] *= sign;
            this->Dzx[idx] *= sign;
            this->Dzz[idx] *= sign;
        }
    }
}

inline float calculate_norm(const float x, const float y, const float z) 
{
    return sqrt(x*x + y*y + z*z);
}

void FFTOcean::form_xyz_array(void) {
    bool calculate_choppy = this->choppy_lambda > 1e-6f ? true : false;
    const float lambda = this->choppy_lambda;
    for(size_t i = 0; i < this->N; i++) {
        for(size_t j = 0; j < this->N; j++) {
            size_t idx = i * N + j;
            float x = ((float)i / (float)(N-1)) * Length - Length / 2.0f;
            float z = ((float)j / (float)(N-1)) * Length - Length / 2.0f;
            float y = this->ht[idx].real();

            this->xyz_coord[idx*3+0] = x;
            this->xyz_coord[idx*3+1] = y;
            this->xyz_coord[idx*3+2] = z;

            float dx = this->Gx[idx].real();
            float dz = this->Gz[idx].real();
            if (calculate_choppy) {
                // calculate horizontal displacement
                float choppy_x = this->Dx[idx].real();
                float choppy_z = this->Dz[idx].real();
                this->xyz_coord[idx*3+0] += this->choppy_lambda * choppy_x;
                this->xyz_coord[idx*3+2] += this->choppy_lambda * choppy_z;

                // calculate grad
                float nx = dz * lambda * Dzx[idx].real() - (1 + lambda * Dzz[idx].real()) * dx;
                float ny = (1 + lambda * Dzz[idx].real()) * (1 + lambda * Dxx[idx].real()) - lambda * lambda * Dxz[idx].real() * Dzx[idx].real();
                float nz = dx * lambda * Dxz[idx].real() - dz * (1 + lambda * Dxx[idx].real());
                float norm = calculate_norm(nx, ny, nz);
                this->gxyz_coord[idx*3+0] = nx / norm;
                this->gxyz_coord[idx*3+1] = ny / norm;
                this->gxyz_coord[idx*3+2] = nz / norm;
            } else {
                // calculate only grad.
                float grad_norm = sqrt(dx*dx + 1.f*1.f + dz*dz);
                this->gxyz_coord[idx*3+0] = -dx / grad_norm;
                this->gxyz_coord[idx*3+1] = 1.f / grad_norm;
                this->gxyz_coord[idx*3+2] = -dz / grad_norm;
            }
        }
    }
}

void FFTOcean::Update(float t)
{
    this->t = t;
    this->calculate_ht(t);
    this->calculate_grad();
    this->iFFT();
    this->processSign();
    this->form_xyz_array();
}

}   // namespace FFTOcean