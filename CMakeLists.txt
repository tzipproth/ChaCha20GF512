cmake_minimum_required(VERSION 3.16)

project(ChaCha20GF512 VERSION 0.1 LANGUAGES CXX)

option(
    CHACHA20GF512_FORCE_PORTABLE
    "Disable PCLMUL acceleration and force the portable GF(2^64) multiplication path"
    OFF
)

add_executable(ChaCha20GF512
    main.cpp
)

target_compile_features(ChaCha20GF512 PRIVATE cxx_std_17)

find_package(Threads REQUIRED)
target_link_libraries(ChaCha20GF512 PRIVATE Threads::Threads)

if(WIN32)
    # Required by ChaCha20GF512Seed.h for BCryptGenRandom.
    target_link_libraries(ChaCha20GF512 PRIVATE bcrypt)
endif()

if(CHACHA20GF512_FORCE_PORTABLE)
    target_compile_definitions(ChaCha20GF512 PRIVATE CHACHA20GF512_FORCE_PORTABLE=1)
endif()

if(MSVC)
    target_compile_options(ChaCha20GF512 PRIVATE /O2)
else()
    target_compile_options(ChaCha20GF512 PRIVATE -O3)
endif()
