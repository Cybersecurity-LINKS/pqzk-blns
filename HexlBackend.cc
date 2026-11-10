#include "HexlBackend.h"
#include <limits>
#include <stdexcept>
#include <utility>

using NTL::ZZ;
using NTL::zz_p;
using NTL::zz_pX;

namespace
{
    std::uint64_t centered_residue(std::uint64_t coefficient, std::uint64_t q, std::uint64_t p) 
    {
        if (coefficient <= q / 2) return coefficient % p;
        const std::uint64_t magnitude = q - coefficient;
        const std::uint64_t residue = magnitude % p;
        return residue == 0 ? 0 : p - residue;
    }
} // namespace

HexlBackend::HexlBackend(
    std::size_t degree,
    std::uint64_t logical_modulus,
    std::vector<std::uint64_t> crt_moduli
)
    : degree_(degree),
        q_(logical_modulus),
        q_zz_(
            NTL::conv<ZZ>(
                static_cast<long>(logical_modulus)
            )
        ),
        moduli_(std::move(crt_moduli)),
        crt_product_(ZZ(1))
{
    if (degree_ == 0)
    {
        throw std::invalid_argument(
            "polynomial degree must be non-zero"
        );
    }

    if ((degree_ & (degree_ - 1)) != 0)
    {
        throw std::invalid_argument(
            "polynomial degree must be a power of two"
        );
    }

    if (q_ == 0)
    {
        throw std::invalid_argument(
            "logical modulus must be non-zero"
        );
    }

    if (moduli_.empty())
    {
        throw std::invalid_argument(
            "at least one CRT modulus is required"
        );
    }

    configuration_.degree =
    degree_;

    configuration_.logical_modulus =
        q_;

    configuration_.crt_moduli =
        moduli_;

    ntt_.resize(moduli_.size());

    crt_partial_products_.resize(
        moduli_.size()
    );

    crt_partial_inverses_.resize(
        moduli_.size()
    );

    for (std::size_t r = 0;
         r < moduli_.size();
         ++r)
    {
        const std::uint64_t modulus =
            moduli_[r];

        if ((modulus - 1) % (2 * degree_) != 0)
        {
            throw std::invalid_argument(
                "CRT modulus is incompatible with "
                "the requested negacyclic NTT degree"
            );
        }

        ntt_[r] =
            std::make_unique<intel::hexl::NTT>(
                degree_,
                modulus
            );

        crt_product_ *=
            NTL::conv<ZZ>(
                static_cast<long>(modulus)
            );
    }

    crt_half_product_ =
        crt_product_ / 2;

    for (std::size_t r = 0;
         r < moduli_.size();
         ++r)
    {
        const ZZ modulus =
            NTL::conv<ZZ>(
                static_cast<long>(moduli_[r])
            );

        crt_partial_products_[r] =
            crt_product_ / modulus;

        crt_partial_inverses_[r] =
            NTL::InvMod(
                crt_partial_products_[r]
                    % modulus,
                modulus
            );
    }
}

std::uint64_t HexlBackend::logical_modulus() const noexcept { return q_; }

const HexlBackend::Configuration&
HexlBackend::configuration() const noexcept
{
    return configuration_;
}

std::uint64_t HexlBackend::scalar_to_crt(
    const NTL::zz_p& scalar,
    std::size_t crt_index
) const
{
    if (crt_index >= moduli_.size())
    {
        throw std::out_of_range(
            "CRT scalar index is out of range"
        );
    }

    const std::uint64_t representative =
        static_cast<std::uint64_t>(
            NTL::rep(
                scalar
            )
        );

    return centered_residue(
        representative,
        q_,
        moduli_[crt_index]
    );
}

void HexlBackend::
pointwise_multiply_accumulate_with_scratch(
    CrtPoly& output_ntt,
    const CrtPoly& left_ntt,
    const CrtPoly& right_ntt,
    Poly& scratch
) const
{
    if (output_ntt.empty())
    {
        output_ntt.assign(
            moduli_.size(),
            Poly(degree_, 0)
        );
    }

    if (scratch.size() != degree_)
    {
        scratch.resize(
            degree_
        );
    }

    for (std::size_t crt_index = 0;
         crt_index < moduli_.size();
         ++crt_index)
    {
        intel::hexl::EltwiseMultMod(
            scratch.data(),
            left_ntt[crt_index].data(),
            right_ntt[crt_index].data(),
            degree_,
            moduli_[crt_index],
            1
        );

        intel::hexl::EltwiseAddMod(
            output_ntt[crt_index].data(),
            output_ntt[crt_index].data(),
            scratch.data(),
            degree_,
            moduli_[crt_index]
        );
    }
}

void HexlBackend::check_ntl_modulus() const 
{
    if (static_cast<std::uint64_t>(zz_p::modulus()) != q_)
        throw std::logic_error("NTL modulus does not match HexlBackend modulus");
}

void HexlBackend::validate_vector(
    const PreparedVector& vector
) const
{
    if (vector.configuration != configuration())
    {
        throw std::invalid_argument(
            "prepared vector is incompatible with this HexlBackend configuration"
        );
    }
}

void HexlBackend::validate_matrix(const PreparedMatrix& matrix) const
{
    if (matrix.configuration != configuration())
    {
        throw std::invalid_argument(
            "prepared matrix is incompatible with this HexlBackend configuration "
        );
    }

    if (matrix.entries_ntt.size()
        != matrix.rows * matrix.cols)
    {
        throw std::invalid_argument(
            "invalid prepared matrix storage"
        );
    }
}

ZZ HexlBackend::reconstruct_centered(
    const CrtPoly& input,
    std::size_t coefficient_index
) const
{
    if (input.size() != moduli_.size())
    {
        throw std::invalid_argument(
            "invalid number of CRT limbs"
        );
    }

    ZZ x(0);

    for (std::size_t r = 0;
         r < moduli_.size();
         ++r)
    {
        if (input[r].size() != degree_)
        {
            throw std::invalid_argument(
                "invalid CRT polynomial degree"
            );
        }

        x += NTL::conv<ZZ>(
                 static_cast<long>(
                     input[r][coefficient_index]
                 )
             )
           * crt_partial_products_[r]
           * crt_partial_inverses_[r];
    }

    x %= crt_product_;

    if (x < 0)
    {
        x += crt_product_;
    }

    if (x > crt_half_product_)
    {
        x -= crt_product_;
    }

    return x;
}

void HexlBackend::to_crt(CrtPoly& output, const zz_pX& input) const 
{
    check_ntl_modulus();
    output.assign(
        moduli_.size(),
        Poly(degree_, 0)
    );
    const long degree = NTL::deg(input);
    for (std::size_t j = 0; j < degree_; ++j)
    {
        std::uint64_t coefficient = 0;
        if (static_cast<long>(j) <= degree)
            coefficient = static_cast<std::uint64_t>(
                NTL::rep(NTL::coeff(input, static_cast<long>(j))));
        for (std::size_t r = 0; r < moduli_.size(); ++r)
            output[r][j] = centered_residue(coefficient, q_, moduli_[r]);
    }
}

void HexlBackend::from_crt(zz_pX& output, const CrtPoly& input) const 
{
    check_ntl_modulus();
    output.SetLength(static_cast<long>(degree_));
    for (std::size_t j = 0; j < degree_; ++j)
    {
        ZZ x = reconstruct_centered(input, j );
        x %= q_zz_;
        if (x < 0) x += q_zz_;
        output[static_cast<long>(j)] = NTL::conv<zz_p>(x);
    }
    output.normalize();
}

void HexlBackend::forward(CrtPoly& output_ntt, const CrtPoly& input_coeff) const 
{
    output_ntt.assign(
        moduli_.size(),
        Poly(degree_, 0)
    );
    for (std::size_t r = 0; r < moduli_.size(); ++r)
    {
        ntt_[r]->ComputeForward(
            output_ntt[r].data(),
            input_coeff[r].data(),
            1,
            1
        );
    }
}

void HexlBackend::inverse(CrtPoly& output_coeff, const CrtPoly& input_ntt) const 
{
    output_coeff.assign(
        moduli_.size(),
        Poly(degree_, 0)
    );
    for (std::size_t r = 0; r < moduli_.size(); ++r)
    {
        ntt_[r]->ComputeInverse(
            output_coeff[r].data(),
            input_ntt[r].data(),
            1,
            1
        );
    }
}

void HexlBackend::pointwise_multiply(CrtPoly& output_ntt, const CrtPoly& left_ntt, const CrtPoly& right_ntt) const 
{
    output_ntt.assign(
        moduli_.size(),
        Poly(degree_, 0)
    );
    for (std::size_t r = 0; r < moduli_.size(); ++r)
    {
        intel::hexl::EltwiseMultMod(
            output_ntt[r].data(),
            left_ntt[r].data(),
            right_ntt[r].data(),
            degree_,
            moduli_[r],
            1
        );
    }
}

void HexlBackend::pointwise_multiply_accumulate(
    CrtPoly& output_ntt,
    const CrtPoly& left_ntt,
    const CrtPoly& right_ntt
) const
{
    if (output_ntt.empty())
    {
        output_ntt.assign(
            moduli_.size(),
            Poly(degree_, 0)
        );
    }

    Poly product(
        degree_,
        0
    );

    for (std::size_t r = 0;
         r < moduli_.size();
         ++r)
    {
        intel::hexl::EltwiseMultMod(
            product.data(),
            left_ntt[r].data(),
            right_ntt[r].data(),
            degree_,
            moduli_[r],
            1
        );

        intel::hexl::EltwiseAddMod(
            output_ntt[r].data(),
            output_ntt[r].data(),
            product.data(),
            degree_,
            moduli_[r]
        );
    }
}

void HexlBackend::finalize_accumulator(
    zz_pX& output, const CrtPoly& accumulator_ntt) const
{
    CrtPoly accumulator_coeff{};
    inverse(accumulator_coeff, accumulator_ntt);
    from_crt(output, accumulator_coeff);
}

void HexlBackend::multiply(zz_pX& output,
                           const zz_pX& left,
                           const zz_pX& right) const 
{
    CrtPoly left_crt{}, right_crt{}, left_ntt{}, right_ntt{};
    CrtPoly product_ntt{}, product_crt{};
    to_crt(left_crt, left);
    to_crt(right_crt, right);
    forward(left_ntt, left_crt);
    forward(right_ntt, right_crt);
    pointwise_multiply(product_ntt, left_ntt, right_ntt);
    inverse(product_crt, product_ntt);
    from_crt(output, product_crt);
}

void HexlBackend::prepare_vector(PreparedVector& output,
                                 const NTL::Vec<zz_pX>& input) const 
{
    check_ntl_modulus();
    output.configuration = configuration();
    output.entries_ntt.resize(static_cast<std::size_t>(input.length()));
    for (long i = 0; i < input.length(); ++i) {
        CrtPoly coeff{};
        to_crt(coeff, input[i]);
        forward(output.entries_ntt[static_cast<std::size_t>(i)], coeff);
    }
}

void HexlBackend::prepare_matrix(PreparedMatrix& output,
                                 const NTL::Mat<zz_pX>& input) const 
{
    check_ntl_modulus();
    output.configuration = configuration();
    output.rows = static_cast<std::size_t>(input.NumRows());
    output.cols = static_cast<std::size_t>(input.NumCols());
    output.entries_ntt.resize(output.rows * output.cols);
    for (std::size_t row = 0; row < output.rows; ++row) {
        for (std::size_t col = 0; col < output.cols; ++col) {
            CrtPoly coeff{};
            to_crt(coeff, input[static_cast<long>(row)][static_cast<long>(col)]);
            forward(output.entries_ntt[row * output.cols + col], coeff);
        }
    }
}

void HexlBackend::prepare_coefficient_block_matrix(
    PreparedMatrix& output,
    const NTL::Mat<NTL::zz_p>& input,
    std::size_t polynomial_columns
) const
{
    prepare_coefficient_block_matrix(
        output,
        input,
        0,
        polynomial_columns
    );
}

void HexlBackend::prepare_coefficient_block_matrix(
    PreparedMatrix& output,
    const NTL::Mat<NTL::zz_p>& input,
    std::size_t input_column_offset,
    std::size_t polynomial_columns
) const
{
    check_ntl_modulus();

    const std::size_t rows =
        static_cast<std::size_t>(
            input.NumRows()
        );

    const std::size_t input_columns =
        static_cast<std::size_t>(
            input.NumCols()
        );

    if (polynomial_columns == 0)
    {
        throw std::invalid_argument(
            "polynomial column count must be non-zero"
        );
    }

    if (
        polynomial_columns
        > std::numeric_limits<std::size_t>::max()
          / degree_
    )
    {
        throw std::overflow_error(
            "coefficient-block matrix size overflow"
        );
    }

    const std::size_t coefficient_count =
        polynomial_columns * degree_;

    if (
        input_column_offset > input_columns
        || coefficient_count
           > input_columns - input_column_offset
    )
    {
        throw std::invalid_argument(
            "coefficient-block matrix range is out of bounds"
        );
    }

    if (
        rows != 0
        && polynomial_columns
           > std::numeric_limits<std::size_t>::max()
             / rows
    )
    {
        throw std::overflow_error(
            "prepared matrix entry count overflow"
        );
    }

    output.configuration =
        configuration();

    output.rows =
        rows;

    output.cols =
        polynomial_columns;

    output.entries_ntt.resize(
        rows * polynomial_columns
    );

    // Scratch allocated once and completely overwritten
    // for every coefficient block.
    CrtPoly coefficient_crt(
        moduli_.size(),
        Poly(degree_)
    );

    for (std::size_t row = 0;
         row < rows;
         ++row)
    {
        const NTL::Vec<NTL::zz_p>& input_row =
            input[
                static_cast<long>(row)
            ];

        for (std::size_t block = 0;
             block < polynomial_columns;
             ++block)
        {
            const std::size_t block_offset =
                input_column_offset
                + block * degree_;

            for (std::size_t coefficient = 0;
                 coefficient < degree_;
                 ++coefficient)
            {
                const std::uint64_t representative =
                    static_cast<std::uint64_t>(
                        NTL::rep(
                            input_row[
                                static_cast<long>(
                                    block_offset
                                    + coefficient
                                )
                            ]
                        )
                    );

                for (std::size_t crt_index = 0;
                     crt_index < moduli_.size();
                     ++crt_index)
                {
                    coefficient_crt[
                        crt_index
                    ][
                        coefficient
                    ] = centered_residue(
                        representative,
                        q_,
                        moduli_[crt_index]
                    );
                }
            }

            forward(
                output.entries_ntt[
                    row * polynomial_columns
                    + block
                ],
                coefficient_crt
            );
        }
    }
}

void HexlBackend::automorphism_ntt_inplace(
    CrtPoly& value_ntt
) const
{
    if (value_ntt.size() != moduli_.size())
    {
        throw std::invalid_argument(
            "invalid CRT limb count in automorphism"
        );
    }

    for (std::size_t r = 0;
         r < moduli_.size();
         ++r)
    {
        if (value_ntt[r].size() != degree_)
        {
            throw std::invalid_argument(
                "invalid polynomial degree in automorphism"
            );
        }

        for (std::size_t j = 0;
             j < degree_ / 2;
             ++j)
        {
            std::swap(
                value_ntt[r][j],
                value_ntt[r][degree_ - 1 - j]
            );
        }
    }
}

void HexlBackend::automorphism_ntt(
    CrtPoly& output_ntt,
    const CrtPoly& input_ntt
) const
{
    if (&output_ntt == &input_ntt)
    {
        automorphism_ntt_inplace(
            output_ntt
        );

        return;
    }

    if (input_ntt.size() != moduli_.size())
    {
        throw std::invalid_argument(
            "invalid CRT limb count in automorphism"
        );
    }

    output_ntt.assign(
        moduli_.size(),
        Poly(degree_, 0)
    );

    for (std::size_t r = 0;
         r < moduli_.size();
         ++r)
    {
        if (input_ntt[r].size() != degree_)
        {
            throw std::invalid_argument(
                "invalid polynomial degree in automorphism"
            );
        }

        for (std::size_t j = 0;
             j < degree_;
             ++j)
        {
            output_ntt[r][j] =
                input_ntt[r][degree_ - 1 - j];
        }
    }
}

void HexlBackend::automorphism_inplace(
    PreparedVector& value
) const
{
    validate_vector(value);

    for (CrtPoly& entry : value.entries_ntt)
    {
        automorphism_ntt_inplace(
            entry
        );
    }
}

void HexlBackend::automorphism(
    PreparedVector& output,
    const PreparedVector& input
) const
{
    validate_vector(input);

    if (&output == &input)
    {
        automorphism_inplace(
            output
        );

        return;
    }

    output.configuration =
        configuration();

    output.entries_ntt.resize(
        input.entries_ntt.size()
    );

    for (std::size_t i = 0;
         i < input.entries_ntt.size();
         ++i)
    {
        automorphism_ntt(
            output.entries_ntt[i],
            input.entries_ntt[i]
        );
    }
}

void HexlBackend::automorphism_inplace(
    PreparedMatrix& value
) const
{
    validate_matrix(value);

    for (CrtPoly& entry : value.entries_ntt)
    {
        automorphism_ntt_inplace(
            entry
        );
    }
}

void HexlBackend::automorphism(
    PreparedMatrix& output,
    const PreparedMatrix& input
) const
{
    validate_matrix(input);

    if (&output == &input)
    {
        automorphism_inplace(
            output
        );

        return;
    }

    output.configuration =
        configuration();

    output.rows =
        input.rows;

    output.cols =
        input.cols;

    output.entries_ntt.resize(
        input.entries_ntt.size()
    );

    for (std::size_t i = 0;
         i < input.entries_ntt.size();
         ++i)
    {
        automorphism_ntt(
            output.entries_ntt[i],
            input.entries_ntt[i]
        );
    }
}

void HexlBackend::matrix_vector_multiply_prepared(
    NTL::Vec<zz_pX>& output,
    const PreparedMatrix& matrix,
    const PreparedVector& vector) const 
{
    validate_matrix(matrix);
    validate_vector(vector);

    if (vector.size() != matrix.cols)
        throw std::invalid_argument("matrix-vector dimension mismatch");

    output.SetLength(static_cast<long>(matrix.rows));
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        CrtPoly accumulator_ntt(
            moduli_.size(),
            Poly(degree_, 0)
        );
        for (std::size_t col = 0; col < matrix.cols; ++col)
            pointwise_multiply_accumulate(accumulator_ntt, matrix.at(row, col),
                                          vector.entries_ntt[col]);
        finalize_accumulator(output[static_cast<long>(row)], accumulator_ntt);
    }
}

void HexlBackend::matrix_vector_multiply_to_prepared(
    PreparedVector& output_ntt,
    const PreparedMatrix& matrix_ntt,
    const PreparedVector& vector_ntt
) const
{
    validate_matrix(
        matrix_ntt
    );

    validate_vector(
        vector_ntt
    );

    if (
        matrix_ntt.cols
        != vector_ntt.size()
    )
    {
        throw std::invalid_argument(
            "matrix-vector dimension mismatch"
        );
    }

    output_ntt.configuration =
        configuration();

    output_ntt.entries_ntt.resize(
        matrix_ntt.rows
    );

    Poly scratch(
        degree_
    );

    for (std::size_t row = 0;
         row < matrix_ntt.rows;
         ++row)
    {
        CrtPoly accumulator_ntt(
            moduli_.size(),
            Poly(degree_, 0)
        );

        for (std::size_t column = 0;
             column < matrix_ntt.cols;
             ++column)
        {
            pointwise_multiply_accumulate_with_scratch(
                accumulator_ntt,
                matrix_ntt.at(
                    row,
                    column
                ),
                vector_ntt.entries_ntt[
                    column
                ],
                scratch
            );
        }

        output_ntt.entries_ntt[row] =
            std::move(
                accumulator_ntt
            );
    }
}
void HexlBackend::combine_three_vectors_with_scalar_matrix(
    PreparedVector& output_ntt,
    const NTL::Mat<NTL::zz_p>& scalars,
    std::size_t scalar_column_offset,
    const PreparedVector& positive,
    const PreparedVector& negative_first,
    const PreparedVector& negative_second
) const
{
    check_ntl_modulus();

    validate_vector(
        positive
    );

    validate_vector(
        negative_first
    );

    validate_vector(
        negative_second
    );

    const std::size_t inner_dimension =
        positive.size();

    if (
        negative_first.size() != inner_dimension
        || negative_second.size() != inner_dimension
    )
    {
        throw std::invalid_argument(
            "prepared vector length mismatch"
        );
    }

    const std::size_t output_rows =
        static_cast<std::size_t>(
            scalars.NumRows()
        );

    const std::size_t scalar_columns =
        static_cast<std::size_t>(
            scalars.NumCols()
        );

    if (
        scalar_column_offset > scalar_columns
        || inner_dimension
           > scalar_columns - scalar_column_offset
    )
    {
        throw std::invalid_argument(
            "scalar matrix range is out of bounds"
        );
    }

    output_ntt.configuration =
        configuration();

    output_ntt.entries_ntt.resize(
        output_rows
    );

    for (std::size_t row = 0;
         row < output_rows;
         ++row)
    {
        CrtPoly accumulator(
            moduli_.size(),
            Poly(degree_, 0)
        );

        for (std::size_t column = 0;
             column < inner_dimension;
             ++column)
        {
            const NTL::zz_p& scalar =
                scalars[
                    static_cast<long>(row)
                ][
                    static_cast<long>(
                        scalar_column_offset
                        + column
                    )
                ];

            for (std::size_t crt_index = 0;
                 crt_index < moduli_.size();
                 ++crt_index)
            {
                const std::uint64_t modulus =
                    moduli_[crt_index];

                const std::uint64_t scalar_crt =
                    scalar_to_crt(
                        scalar,
                        crt_index
                    );

                std::uint64_t* accumulator_data =
                    accumulator[
                        crt_index
                    ].data();

                const std::uint64_t* positive_data =
                    positive.entries_ntt[
                        column
                    ][
                        crt_index
                    ].data();

                const std::uint64_t* negative_first_data =
                    negative_first.entries_ntt[
                        column
                    ][
                        crt_index
                    ].data();

                const std::uint64_t* negative_second_data =
                    negative_second.entries_ntt[
                        column
                    ][
                        crt_index
                    ].data();

                for (std::size_t slot = 0;
                     slot < degree_;
                     ++slot)
                {
                    // combined =
                    //     positive
                    //     - negative_first
                    //     - negative_second
                    // mod p_r.
                    std::uint64_t combined =
                        positive_data[slot];

                    combined =
                        combined
                            >= negative_first_data[slot]
                        ? combined
                          - negative_first_data[slot]
                        : modulus
                          - (
                              negative_first_data[slot]
                              - combined
                          );

                    combined =
                        combined
                            >= negative_second_data[slot]
                        ? combined
                          - negative_second_data[slot]
                        : modulus
                          - (
                              negative_second_data[slot]
                              - combined
                          );

                    const std::uint64_t product =
                        static_cast<std::uint64_t>(
                            (
                                static_cast<
                                    unsigned __int128
                                >(
                                    scalar_crt
                                )
                                * combined
                            )
                            % modulus
                        );

                    const std::uint64_t current =
                        accumulator_data[slot];

                    accumulator_data[slot] =
                        current >= modulus - product
                        ? current
                          - (modulus - product)
                        : current + product;
                }
            }
        }

        output_ntt.entries_ntt[row] =
            std::move(
                accumulator
            );
    }
}

void HexlBackend::matrix_vector_multiply(
    NTL::Vec<zz_pX>& output,
    const PreparedMatrix& matrix,
    const NTL::Vec<zz_pX>& vector) const {
    PreparedVector prepared;
    prepare_vector(prepared, vector);
    matrix_vector_multiply_prepared(output, matrix, prepared);
}

void HexlBackend::matrix_vector_multiply_add_prepared(
    NTL::Vec<zz_pX>& output,
    const PreparedMatrix& left_matrix,
    const PreparedVector& left_vector,
    const PreparedMatrix& right_matrix,
    const PreparedVector& right_vector) const 
{
    validate_matrix(left_matrix);
    validate_matrix(right_matrix);

    validate_vector(left_vector);
    validate_vector(right_vector);

    if (left_matrix.rows != right_matrix.rows)
        throw std::invalid_argument("matrix row count mismatch");
    if (left_matrix.cols != left_vector.size() ||
        right_matrix.cols != right_vector.size())
        throw std::invalid_argument("matrix-vector dimension mismatch");

    output.SetLength(static_cast<long>(left_matrix.rows));
    for (std::size_t row = 0; row < left_matrix.rows; ++row) {
        CrtPoly accumulator_ntt(
            moduli_.size(),
            Poly(degree_, 0)
        );
        for (std::size_t col = 0; col < left_matrix.cols; ++col)
            pointwise_multiply_accumulate(accumulator_ntt,
                                          left_matrix.at(row, col),
                                          left_vector.entries_ntt[col]);
        for (std::size_t col = 0; col < right_matrix.cols; ++col)
            pointwise_multiply_accumulate(accumulator_ntt,
                                          right_matrix.at(row, col),
                                          right_vector.entries_ntt[col]);
        finalize_accumulator(output[static_cast<long>(row)], accumulator_ntt);
    }
}

void HexlBackend::matrix_vector_multiply_add(
    NTL::Vec<zz_pX>& output,
    const PreparedMatrix& left_matrix,
    const NTL::Vec<zz_pX>& left_vector,
    const PreparedMatrix& right_matrix,
    const NTL::Vec<zz_pX>& right_vector) const 
{
    PreparedVector left_prepared, right_prepared;
    prepare_vector(left_prepared, left_vector);
    prepare_vector(right_prepared, right_vector);
    matrix_vector_multiply_add_prepared(output, left_matrix, left_prepared,
                                        right_matrix, right_prepared);
}

void HexlBackend::inner_product_prepared(
    zz_pX& output,
    const PreparedVector& left,
    const PreparedVector& right
) const
{
    validate_vector(left);
    validate_vector(right);

    if (left.size() != right.size())
    {
        throw std::invalid_argument(
            "prepared inner-product dimension mismatch"
        );
    }

    CrtPoly accumulator_ntt(
        moduli_.size(),
        Poly(degree_, 0)
    );

    for (std::size_t i = 0;
         i < left.size();
         ++i)
    {
        pointwise_multiply_accumulate(
            accumulator_ntt,
            left.entries_ntt[i],
            right.entries_ntt[i]
        );
    }

    finalize_accumulator(
        output,
        accumulator_ntt
    );
}

void HexlBackend::inner_product_add_prepared(
    zz_pX& output,
    const PreparedVector& left_first,
    const PreparedVector& right_first,
    const PreparedVector& left_second,
    const PreparedVector& right_second
) const
{
    validate_vector(left_first);
    validate_vector(right_first);
    validate_vector(left_second);
    validate_vector(right_second);

    if (left_first.size() != right_first.size())
    {
        throw std::invalid_argument(
            "first prepared inner-product dimension mismatch"
        );
    }

    if (left_second.size() != right_second.size())
    {
        throw std::invalid_argument(
            "second prepared inner-product dimension mismatch"
        );
    }

    CrtPoly accumulator_ntt(
        moduli_.size(),
        Poly(degree_, 0)
    );

    for (std::size_t i = 0;
         i < left_first.size();
         ++i)
    {
        pointwise_multiply_accumulate(
            accumulator_ntt,
            left_first.entries_ntt[i],
            right_first.entries_ntt[i]
        );
    }

    for (std::size_t i = 0;
         i < left_second.size();
         ++i)
    {
        pointwise_multiply_accumulate(
            accumulator_ntt,
            left_second.entries_ntt[i],
            right_second.entries_ntt[i]
        );
    }

    finalize_accumulator(
        output,
        accumulator_ntt
    );
}

void HexlBackend::inner_product_split_prepared(
    zz_pX& output,
    const PreparedVector& left_first,
    const PreparedVector& left_second,
    const PreparedVector& right_combined,
    std::size_t split
) const
{
    validate_vector(left_first);
    validate_vector(left_second);
    validate_vector(right_combined);

    if (left_first.size() != split)
    {
        throw std::invalid_argument(
            "invalid first split size"
        );
    }

    if (right_combined.size()
        != left_first.size()
         + left_second.size())
    {
        throw std::invalid_argument(
            "split inner-product dimension mismatch"
        );
    }

    CrtPoly accumulator_ntt(
        moduli_.size(),
        Poly(degree_, 0)
    );

    for (std::size_t i = 0;
         i < left_first.size();
         ++i)
    {
        pointwise_multiply_accumulate(
            accumulator_ntt,
            left_first.entries_ntt[i],
            right_combined.entries_ntt[i]
        );
    }

    for (std::size_t i = 0;
         i < left_second.size();
         ++i)
    {
        pointwise_multiply_accumulate(
            accumulator_ntt,
            left_second.entries_ntt[i],
            right_combined.entries_ntt[
                split + i
            ]
        );
    }

    finalize_accumulator(
        output,
        accumulator_ntt
    );
}

void HexlBackend::vector_matrix_multiply_prepared(
    PreparedVector& output_ntt,
    const PreparedVector& vector_ntt,
    const PreparedMatrix& matrix_ntt
) const
{
    validate_vector(vector_ntt);
    validate_matrix(matrix_ntt);

    if (vector_ntt.size() != matrix_ntt.rows)
    {
        throw std::invalid_argument(
            "vector-matrix dimension mismatch"
        );
    }

    output_ntt.configuration =
        configuration();

    output_ntt.entries_ntt.resize(
        matrix_ntt.cols
    );

    for (std::size_t col = 0;
         col < matrix_ntt.cols;
         ++col)
    {
        CrtPoly accumulator_ntt(
            moduli_.size(),
            Poly(degree_, 0)
        );

        for (std::size_t row = 0;
             row < matrix_ntt.rows;
             ++row)
        {
            pointwise_multiply_accumulate(
                accumulator_ntt,
                vector_ntt.entries_ntt[row],
                matrix_ntt.at(row, col)
            );
        }

        output_ntt.entries_ntt[col] =
            std::move(accumulator_ntt);
    }
}

void HexlBackend::finalize_vector(
    NTL::Vec<zz_pX>& output,
    const PreparedVector& input_ntt
) const
{
    validate_vector(input_ntt);

    output.SetLength(
        static_cast<long>(
            input_ntt.size()
        )
    );

    for (std::size_t i = 0;
         i < input_ntt.size();
         ++i)
    {
        finalize_accumulator(
            output[static_cast<long>(i)],
            input_ntt.entries_ntt[i]
        );
    }
}