#ifdef __EMSCRIPTEN__
#   include <emscripten.h>
#   include <emscripten/bind.h>
#endif

#include "FFTOcean.hpp"

#ifdef __EMSCRIPTEN__
using namespace FFTOceanNamespace;
EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::class_<FFTOcean>("FFTOcean")
        .constructor<unsigned long, unsigned long, float, float, float, float, float>()
        .function("Update", &FFTOcean::Update)
        .function("get_xyz_ptr", &FFTOcean::get_xyz_ptr)
        .function("get_gxyz_ptr", &FFTOcean::get_gxyz_ptr)
        .function("set_height_scale", &FFTOcean::set_height_scale)
        .function("get_height_scale", &FFTOcean::get_height_scale)
        .function("set_choppy_coefficient", &FFTOcean::set_choppy_coefficient)
        .function("get_choppy_coefficient", &FFTOcean::get_choppy_coefficient)
        .function("get_ht_real_ptr", &FFTOcean::get_ht_real_ptr)
        .function("get_Gx_real_ptr", &FFTOcean::get_Gx_real_ptr, emscripten::allow_raw_pointers())
        .function("get_Gz_real_ptr", &FFTOcean::get_Gz_real_ptr, emscripten::allow_raw_pointers())
        .function("get_Dx_real_ptr", &FFTOcean::get_Dx_real_ptr, emscripten::allow_raw_pointers())
        .function("get_Dz_real_ptr", &FFTOcean::get_Dz_real_ptr, emscripten::allow_raw_pointers())
        ;
}
#else

#endif
