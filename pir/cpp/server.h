//
// Copyright 2020 the authors listed in CONTRIBUTORS.md
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#ifndef PIR_SERVER_H_
#define PIR_SERVER_H_

#include <vector>

#include "absl/status/statusor.h"
#include "pir/cpp/context.h"
#include "pir/cpp/database.h"
#include "pir/cpp/serialization.h"
#include "seal/seal.h"

namespace pir {

using absl::Status;
using absl::StatusOr;
using ::seal::GaloisKeys;
using ::seal::RelinKeys;

class PIRServer {
   public:
    /**
     * Creates and returns a new server instance, holding a database.
     * @param[in] db PIRDatabase to load
     * @param[in] params PIR Paramerters
     * @returns InvalidArgument if the database encoding fails
     **/
    static StatusOr<std::unique_ptr<PIRServer>> Create(
        std::shared_ptr<PIRDatabase> database,
        shared_ptr<PIRParameters> params);

    /**
     * Handles a client request.
     * @param[in] request The PIR Payload
     * @returns InvalidArgument if the deserialization or encrypted operations
     *fail
     **/
    StatusOr<Response> ProcessRequest(const Request& request) const;

    PIRServer() = delete;

    /**
     * Helper function to do the substitution operation on a ciphertext. If the
     * ciphertext is the encryption of polynomial p(x), then given power k, the
     * result will be p(x**k).
     * @param ct Ciphertext to operate, modification done in place.
     * @param[in] power The power k to raise x to in the plaintext polynomial.
     * @param[in] gal_keys Galois keys used for automorphism. Must be generated
     * by whoever encrypted the ciphertext using keygen, and must include the
     * power being asked for.
     */
    Status substitute_power_x_inplace(seal::Ciphertext& ct, std::uint32_t power,
                                      const seal::GaloisKeys& gal_keys) const;

    /**
     * Helper function to multiply a ciphertext by a given power of 1/x. As a
     * result plaintext is also multiplied by the same power of 1/x. For
     * example, if the ciphertext is the encryption of p(x) = 99x^5, and if the
     * given k is 3, then this results in p(x) * 1/x^3 = 99x^2.
     * @param[in] encrypted Ciphertext to take as input.
     * @param[in] k Power of 1/x to multiply.
     * @param[out] destination Output ciphertext after multiplying by power of
     * x.
     */
    void multiply_inverse_power_of_x(const seal::Ciphertext& encrypted,
                                     uint32_t k,
                                     seal::Ciphertext& destination) const;

    /**
     * Performs an oblivious expansion on an input ciphertext to a vector of
     * ciphertexts. If the input ciphertext is the encryption of a plaintext
     * polynomial of the form a0 + a1*x + a2*x^2 + ... + an*x^n, where n is the
     * num_items below, then the output is a series of ciphertexts, where each
     * is the encryption of each term as the constant coefficient. In other
     * words, the output will be a vector of: [enc(a0), enc(a1), ..., enc(an)],
     * where enc(ai) is the encryption of a polynomial that only has ai in the
     * constant coefficient (all other terms are zero).
     *
     * The most common example of this is to expand a selection vector in PIR to
     * individual ciphertexts. The selection vector is just 0 in all slots
     * except for 1 in the desired slot. If there are 4 items, and the third
     * item is the one desired, the selection vector is [0, 0, 1, 0]. The client
     * represents this as a single polynomial of the form 1*x^2. The server uses
     * this expansion function to expand this to 4 ciphertexts: [enc(0), enc(0),
     * enc(1), enc(0)], which it then uses to produce the response.
     *
     * NB: Due to an optimization in PIR, this is left as (expanded_vector) * m,
     * where m is the smallest power of 2 greater than num_items. It is assumed
     * that the plaintext modulus will be changed to make this irrelevant.
     *
     * @param[in] ct The input ciphertext to expand.
     * @param[in] num_items The number of items to extract.
     * @param[in] gal_keys Galois keys supplied by the client.
     * @returns A vector of ciphertexts that are the expansion as described
     * above.
     */
    StatusOr<std::vector<seal::Ciphertext>> oblivious_expansion(
        const seal::Ciphertext& ct, const size_t num_items,
        const seal::GaloisKeys& gal_keys) const;

    /**
     * Extension of oblivious_expansion to multiple ciphertexts. This allows
     * selection vectors that are larger then poly_modulus_degree to be used.
     * The output of the expansion of each ciphertext is concatenated to form
     * the results of this function. Each ciphertext that isn't the last in the
     * vector are assumed to contain exactly poly_modulus_degree items to be
     * expanded.
     *
     * @param[in] cts List of ciphertexts to use as input to the expansion
     * @param[in] total_items Total number of ciphertexts after expansion
     * @param[in] gal_keys Galois keys supplied by the client
     * @returns A vector of ciphertexts that are the expansion of all of the
     * input ciphertexts concatenated.
     */
    StatusOr<std::vector<seal::Ciphertext>> oblivious_expansion(
        const std::vector<seal::Ciphertext>& cts, const size_t total_items,
        const seal::GaloisKeys& gal_keys) const;

    // Just for testing: get the context
    PIRContext* Context() { return context_.get(); }

   private:
    PIRServer(std::unique_ptr<PIRContext> /*sealctx*/,
              std::shared_ptr<PIRDatabase> /*db*/);

    Status processQuery(const Ciphertexts& query, const GaloisKeys& galois_keys,
                        const optional<RelinKeys>& relin_keys,
                        const size_t& dim_sum, Ciphertexts* output) const;

    std::unique_ptr<PIRContext> context_;
    std::shared_ptr<PIRDatabase> db_;
};

}  // namespace pir

#endif  // PIR_SERVER_H_
