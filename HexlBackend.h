#ifndef HEXL_BACKEND_H
#define HEXL_BACKEND_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <NTL/ZZ.h>
#include <NTL/lzz_pX.h>
#include <hexl/hexl.hpp>
#include "params.h"

class HexlBackend {
public:

    using Poly = std::vector<std::uint64_t>;
    using CrtPoly = std::vector<Poly>;

    struct Configuration
    {
        std::size_t degree = 0;
        std::uint64_t logical_modulus = 0;
        std::vector<std::uint64_t> crt_moduli;

        bool operator==(
            const Configuration& other
        ) const
        {
            return degree == other.degree
                && logical_modulus == other.logical_modulus
                && crt_moduli == other.crt_moduli;
        }

        bool operator!=(
            const Configuration& other
        ) const
        {
            return !(*this == other);
        }
    };

    struct PreparedVector
    {
        Configuration configuration;

        std::vector<CrtPoly> entries_ntt;

        std::size_t size() const noexcept
        {
            return entries_ntt.size();
        }
    };

    struct PreparedMatrix
    {
        Configuration configuration;

        std::size_t rows = 0;
        std::size_t cols = 0;

        std::vector<CrtPoly> entries_ntt;

        const CrtPoly& at(
            std::size_t row,
            std::size_t col
        ) const
        {
            return entries_ntt.at(
                row * cols + col
            );
        }
    };

    HexlBackend(
        std::size_t degree,
        std::uint64_t logical_modulus,
        std::vector<std::uint64_t> crt_moduli
    );
    

    std::uint64_t logical_modulus() const noexcept;
    const Configuration& configuration() const noexcept;

    void to_crt(CrtPoly& output, const NTL::zz_pX& input) const;
    void from_crt(NTL::zz_pX& output, const CrtPoly& input) const;
    void forward(CrtPoly& output_ntt, const CrtPoly& input_coeff) const;
    void inverse(CrtPoly& output_coeff, const CrtPoly& input_ntt) const;

    void pointwise_multiply(
        CrtPoly& output_ntt,
        const CrtPoly& left_ntt,
        const CrtPoly& right_ntt) const;

    void pointwise_multiply_accumulate(
        CrtPoly& output_ntt,
        const CrtPoly& left_ntt,
        const CrtPoly& right_ntt) const;

    void multiply(
        NTL::zz_pX& output,
        const NTL::zz_pX& left,
        const NTL::zz_pX& right) const;

    void prepare_vector(
        PreparedVector& output,
        const NTL::Vec<NTL::zz_pX>& input) const;

    void prepare_matrix(
        PreparedMatrix& output,
        const NTL::Mat<NTL::zz_pX>& input
    ) const;

    void prepare_coefficient_block_matrix(
        PreparedMatrix& output,
        const NTL::Mat<NTL::zz_p>& input,
        std::size_t polynomial_columns
    ) const;

    void prepare_coefficient_block_matrix(
        PreparedMatrix& output,
        const NTL::Mat<NTL::zz_p>& input,
        std::size_t input_column_offset,
        std::size_t polynomial_columns
    ) const;

    void automorphism_ntt_inplace(
        CrtPoly& value_ntt) const;

    void automorphism_ntt(
        CrtPoly& output_ntt,
        const CrtPoly& input_ntt) const;

    void automorphism_inplace(
        PreparedVector& value) const;

    void automorphism(
        PreparedVector& output,
        const PreparedVector& input) const;

    void automorphism_inplace(
        PreparedMatrix& value) const;

    void automorphism(
        PreparedMatrix& output,
        const PreparedMatrix& input) const;
    
    void inner_product_prepared(
        NTL::zz_pX& output,
        const PreparedVector& left,
        const PreparedVector& right) const;

    void inner_product_add_prepared(
        NTL::zz_pX& output,
        const PreparedVector& left_first,
        const PreparedVector& right_first,
        const PreparedVector& left_second,
        const PreparedVector& right_second) const;

    void inner_product_split_prepared(
        NTL::zz_pX& output,
        const PreparedVector& left_first,
        const PreparedVector& left_second,
        const PreparedVector& right_combined,
        std::size_t split) const;

    void matrix_vector_multiply_prepared(
        NTL::Vec<NTL::zz_pX>& output,
        const PreparedMatrix& matrix,
        const PreparedVector& vector) const;

    void matrix_vector_multiply_to_prepared(
        PreparedVector& output_ntt,
        const PreparedMatrix& matrix_ntt,
        const PreparedVector& vector_ntt
    ) const;

    void combine_three_vectors_with_scalar_matrix(
        PreparedVector& output_ntt,
        const NTL::Mat<NTL::zz_p>& scalars,
        std::size_t scalar_column_offset,
        const PreparedVector& positive,
        const PreparedVector& negative_first,
        const PreparedVector& negative_second
    ) const;

    void matrix_vector_multiply(
        NTL::Vec<NTL::zz_pX>& output,
        const PreparedMatrix& matrix,
        const NTL::Vec<NTL::zz_pX>& vector) const;

    void matrix_vector_multiply_add_prepared(
        NTL::Vec<NTL::zz_pX>& output,
        const PreparedMatrix& left_matrix,
        const PreparedVector& left_vector,
        const PreparedMatrix& right_matrix,
        const PreparedVector& right_vector) const;

    void matrix_vector_multiply_add(
        NTL::Vec<NTL::zz_pX>& output,
        const PreparedMatrix& left_matrix,
        const NTL::Vec<NTL::zz_pX>& left_vector,
        const PreparedMatrix& right_matrix,
        const NTL::Vec<NTL::zz_pX>& right_vector) const;
    std::size_t degree() const noexcept
    {
        return degree_;
    }

    void vector_matrix_multiply_prepared(
        PreparedVector& output_ntt,
        const PreparedVector& vector_ntt,
        const PreparedMatrix& matrix_ntt
    ) const;

    void finalize_vector(
        NTL::Vec<NTL::zz_pX>& output,
        const PreparedVector& input_ntt
    ) const;

    std::size_t crt_moduli_count() const noexcept
    {
        return moduli_.size();
    }

private:

    std::size_t degree_;
    std::uint64_t q_;
    NTL::ZZ q_zz_;

    std::vector<std::uint64_t> moduli_;

    Configuration configuration_;

    std::vector<
        std::unique_ptr<intel::hexl::NTT>
    > ntt_;

    NTL::ZZ crt_product_;
    NTL::ZZ crt_half_product_;

    std::vector<NTL::ZZ> crt_partial_products_;
    std::vector<NTL::ZZ> crt_partial_inverses_;

    std::uint64_t scalar_to_crt(
        const NTL::zz_p& scalar,
        std::size_t crt_index
    ) const;

    void pointwise_multiply_accumulate_with_scratch(
        CrtPoly& output_ntt,
        const CrtPoly& left_ntt,
        const CrtPoly& right_ntt,
        Poly& scratch
    ) const;

    void check_ntl_modulus() const;
    void validate_vector(const PreparedVector& vector) const;
    void validate_matrix(const PreparedMatrix& matrix) const;
    void finalize_accumulator(NTL::zz_pX& output,
                              const CrtPoly& accumulator_ntt) const;

    NTL::ZZ reconstruct_centered(
        const CrtPoly& input,
        std::size_t coefficient_index
    ) const;
};

#endif
