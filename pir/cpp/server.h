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

#include "context.h"
#include "database.h"
#include "payload.h"
#include "seal/seal.h"
#include "util/statusor.h"

namespace pir {

using ::private_join_and_compute::StatusOr;

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
      std::shared_ptr<PIRParameters> params);
  /**
   * Creates and returns a new server instance, holding a database.
   * @param[in] db PIRDatabase to load
   * @returns InvalidArgument if the database encoding fails
   **/
  static StatusOr<std::unique_ptr<PIRServer>> Create(
      std::shared_ptr<PIRDatabase> database);

  /**
   * Handles a client request.
   * @param[in] request The PIR Payload
   * @returns InvalidArgument if the deserialization or encrypted operations
   *fail
   **/
  StatusOr<PIRPayload> ProcessRequest(const PIRPayload& request) const;

  /**
   * Returns the database size.
   **/
  std::size_t DBSize() { return context_->DBSize(); }

  PIRServer() = delete;

  /**
   * Helper function to do the substitution operation on a ciphertext. If the
   * ciphertext is the encryption of polynomial p(x), then given power k, the
   * result will be p(x**k).
   * @param ct Ciphertext to operate, modification done in place.
   * @param[in] power The power k to raise x to in the plaintext polynomial.
   * @param[in] gal_keys Galois keys used for automorphism. Must be generated by
   *   whoever encrypted the ciphertext using keygen, and must include the power
   *   being asked for.
   */
  void substitute_power_x_inplace(seal::Ciphertext& ct, std::uint32_t power,
                                  const seal::GaloisKeys& gal_keys);

  /**
   * Helper function to multiply a ciphertext by a given power of x. As a result
   * plaintext is also multiplied by the same power of x. For example, if the
   * ciphertext is the encryption of p(x) = 99x^2, and if the given k is 3, then
   * this results in p(x) * x^3 = 99x^5.
   * @param[in] encrypted Ciphertext to take as input.
   * @param[in] k Power of x to multiply. Can be negative.
   * @param[out] destination Output ciphertext after multiplying by power of x.
   */
  void multiply_power_of_x(const seal::Ciphertext& encrypted, int k,
                           seal::Ciphertext& destination);

  // Just for testing: get the context
  PIRContext* Context() { return context_.get(); }

 private:
  PIRServer(std::unique_ptr<PIRContext> /*sealctx*/,
            std::shared_ptr<PIRDatabase> /*db*/);

  std::unique_ptr<PIRContext> context_;
  std::shared_ptr<PIRDatabase> db_;
};

}  // namespace pir

#endif  // PIR_SERVER_H_
