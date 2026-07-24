#ifndef BLNS_TERNARY_H
#define BLNS_TERNARY_H

#include "params.h"
#include <stdint.h>

//==============================================================================
// TernaryCoeffStructure
//
// Sparse representation of a vector of ternary polynomial blocks.
//
// entries:
//   flat array containing all non-zero coefficients.
//   Each byte packs:
//     - bit 7     : sign   (1 => +1, 0 => -1)
//     - bits 0..6 : coefficient position inside the polynomial block
//
// offsets:
//   offsets[k] .. offsets[k+1]-1 is the range of entries belonging to
//   polynomial block k.
//
// Assumption:
//   d <= 128, so coefficient positions 0..d-1 fit into 7 bits.
//==============================================================================
struct TernaryCoeffStructure
{
    vector<uint8_t> entries;    // [sign | pos]
    vector<uint16_t> offsets;
};

//==============================================================================
// Packed ternary entry helpers
//==============================================================================

static inline uint8_t ternary_pack_entry(const uint8_t pos, const bool is_plus)
{
    assert(pos < 128);
    return static_cast<uint8_t>(pos | (is_plus ? 0x80u : 0u));
}

static inline uint8_t ternary_entry_pos(const uint8_t code)
{
    return static_cast<uint8_t>(code & 0x7Fu);
}

static inline bool ternary_entry_is_plus(const uint8_t code)
{
    return (code & 0x80u) != 0;
}

//==============================================================================
// Sigma index helpers
//==============================================================================

static inline uint8_t ternary_sigma_index_from_raw(
    const uint8_t raw_pos,
    const long d)
{
    assert(d > 0 && d <= 128);
    assert(static_cast<long>(raw_pos) < d);

    return raw_pos == 0
        ? 0
        : static_cast<uint8_t>(d - raw_pos);
}

static inline uint8_t ternary_sigma_index_to_raw(
    const uint8_t sigma_pos,
    const long d)
{
    assert(d > 0 && d <= 128);
    assert(static_cast<long>(sigma_pos) < d);

    return sigma_pos == 0
        ? 0
        : static_cast<uint8_t>(d - sigma_pos);
}

static inline bool ternary_sigma_flips_sign_from_raw(const uint8_t raw_pos)
{
    return raw_pos != 0;
}

static inline bool ternary_sigma_flips_sign_from_sigma(const uint8_t sigma_pos)
{
    return sigma_pos != 0;
}

//==============================================================================
// Append helpers
//==============================================================================

static inline void ternary_append_raw(
    TernaryCoeffStructure& out,
    const uint8_t raw_pos,
    const int val,
    const long d)
{
    assert(d > 0 && d <= 128);
    assert(static_cast<long>(raw_pos) < d);

    if (val == 0)
        return;

    const bool is_plus = val > 0;
    out.entries.push_back(ternary_pack_entry(raw_pos, is_plus));
}

static inline void ternary_append_sigma(
    TernaryCoeffStructure& out,
    const uint8_t raw_pos,
    const int raw_val,
    const long d)
{
    assert(d > 0 && d <= 128);
    assert(static_cast<long>(raw_pos) < d);

    if (raw_val == 0)
        return;

    const uint8_t sigma_pos =
        ternary_sigma_index_from_raw(raw_pos, d);

    const bool raw_plus = raw_val > 0;
    const bool flip =
        ternary_sigma_flips_sign_from_raw(raw_pos);

    const bool sigma_plus = raw_plus ^ flip;

    out.entries.push_back(
        ternary_pack_entry(sigma_pos, sigma_plus)
    );
}

//==============================================================================
// Block helpers
//==============================================================================

static inline void ternary_blocks_init(
    TernaryCoeffStructure& out,
    const size_t n_blocks,
    const size_t reserve_hint)
{
    out.entries.clear();
    out.offsets.resize(n_blocks + 1);
    out.entries.reserve(reserve_hint);
    out.offsets[0] = 0;
}

static inline void ternary_finish_block(
    TernaryCoeffStructure& out,
    const size_t block_idx)
{
    assert(block_idx < out.offsets.size());
    assert(out.entries.size() <= 65535);

    out.offsets[block_idx] =
        static_cast<uint16_t>(out.entries.size());
}

//==============================================================================
// Sigma-packed row dot raw coefficient vector
//
// Used for:
//   R_goth[i] * coeffs_s1
//
// R_goth row is stored sigma-packed.
// coeffs is stored in raw coefficient order.
//==============================================================================
static inline void ternary_sigma_dot_raw_coeffs(
    zz_p& out,
    const TernaryCoeffStructure& row,
    const vec_zz_p& coeffs,
    const ulong n_blocks,
    const long d)
{
    assert(d > 0 && d <= 128);
    assert(row.offsets.size() == n_blocks + 1);
    assert(coeffs.length() >= static_cast<long>(n_blocks * d));

    clear(out);

    const uint8_t* entries = row.entries.data();
    const uint16_t* offs   = row.offsets.data();

    long base = 0;

    for (ulong k = 0; k < n_blocks; k++, base += d)
    {
        const uint16_t begin = offs[k];
        const uint16_t end   = offs[k + 1];

        for (uint16_t u = begin; u < end; u++)
        {
            const uint8_t enc = entries[u];

            const uint8_t sigma_pos = ternary_entry_pos(enc);
            const bool sigma_plus   = ternary_entry_is_plus(enc);

            const uint8_t raw_pos =
                ternary_sigma_index_to_raw(sigma_pos, d);

            const bool raw_plus =
                sigma_plus ^ ternary_sigma_flips_sign_from_sigma(sigma_pos);

            const zz_p& x = coeffs[base + raw_pos];

            if (raw_plus)
                out += x;
            else
                out -= x;
        }
    }
}

//==============================================================================
// Accumulate gamma through a sigma-packed ternary row
//
// Used for constructing A_hat, B_hat, C_hat.
// Since row is already sigma-packed, positions are consumed directly.
//==============================================================================
static inline void ternary_sigma_accumulate_gamma(
    const TernaryCoeffStructure& row,
    const zz_p& gamma_value,
    const vector<zz_p*>& block_reps,
    const ulong n_blocks)
{
    assert(row.offsets.size() == n_blocks + 1);
    assert(block_reps.size() >= n_blocks);

    const zz_p neg_gamma = -gamma_value;
    const zz_p signed_gamma[2] = { neg_gamma, gamma_value };

    const uint8_t* entries = row.entries.data();
    const uint16_t* offs   = row.offsets.data();

    for (ulong k = 0; k < n_blocks; k++)
    {
        zz_p* rep = block_reps[k];

        const uint16_t begin = offs[k];
        const uint16_t end   = offs[k + 1];

        for (uint16_t u = begin; u < end; u++)
        {
            const uint8_t enc = entries[u];

            const uint8_t pos = ternary_entry_pos(enc);
            const bool is_plus = ternary_entry_is_plus(enc);

            rep[pos] += signed_gamma[is_plus ? 1 : 0];
        }
    }
}

static inline void sample_s2_dense_and_packed(
    TernaryCoeffStructure& s2_raw,
    const ulong m2,
    const long d)
{
    assert(d > 0 && d <= 128);
    assert((m2 * d) <= 65535);

    ternary_blocks_init(
        s2_raw,
        m2,
        (2 * m2 * d) / 3 + 128
    );

    for (ulong i = 0; i < m2; i++)
    {

        for (long j = 0; j < d; j++)
        {
            const long r = RandomBnd(3) - 1;

            ternary_append_raw(
                s2_raw,
                static_cast<uint8_t>(j),
                static_cast<int>(r),
                d
            );
        }

        ternary_finish_block(s2_raw, i + 1);
    }
}

#endif