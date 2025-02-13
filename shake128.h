/*
Copyright (c) 2022-2024 IBM

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
 * Based on the public domain implementation in
 * crypto_hash/keccakc512/simple/ from http://bench.cr.yp.to/supercop.html
 * by Ronny Van Keer and the public domain "TweetFips202" implementation
 * from https://twitter.com/tweetfips202 by Gilles Van Assche,
 * Daniel J. Bernstein, and Peter Schwabe
 */

#ifndef SHAKE128_H
#define SHAKE128_H

#include <stddef.h>
#include <stdint.h>


typedef struct {
  uint64_t s[25];
  unsigned int pos;
  int final;
} shake128_state_t;


void _shake128_init (shake128_state_t *state);
void _shake128_absorb (shake128_state_t *state, const uint8_t *in,
                              size_t len);
void _shake128_squeeze (shake128_state_t *state, uint8_t *out,
                               size_t len);

#endif

