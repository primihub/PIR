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
#include "pir/cpp/database.h"

#include <iostream>
#include <memory>
#include <thread>

#include "pir/cpp/ct_reencoder.h"
#include "pir/cpp/status_asserts.h"
#include "pir/cpp/string_encoder.h"
#include "pir/cpp/utils.h"
#include "seal/seal.h"

namespace pir {

using absl::InternalError;
using absl::InvalidArgumentError;
using absl::StatusOr;
using google::protobuf::RepeatedField;
using seal::Ciphertext;
using seal::Evaluator;
using seal::Plaintext;
using std::unique_ptr;
using std::vector;

StatusOr<shared_ptr<PIRDatabase>> PIRDatabase::Create(
    shared_ptr<PIRParameters> params) {
    ASSIGN_OR_RETURN(auto context, PIRContext::Create(params));
    return std::make_shared<PIRDatabase>(std::move(context));
}
StatusOr<shared_ptr<PIRDatabase>> PIRDatabase::Create(
    const vector<std::int64_t>& rawdb, shared_ptr<PIRParameters> params) {
    ASSIGN_OR_RETURN(auto context, PIRContext::Create(params));
    auto pir_db = std::make_shared<PIRDatabase>(std::move(context));
    RETURN_IF_ERROR(pir_db->populate(rawdb));
    return std::move(pir_db);
}

StatusOr<shared_ptr<PIRDatabase>> PIRDatabase::Create(
    const vector<string>& rawdb, shared_ptr<PIRParameters> params) {
    ASSIGN_OR_RETURN(auto context, PIRContext::Create(params));
    auto pir_db = std::make_shared<PIRDatabase>(std::move(context));
    RETURN_IF_ERROR(pir_db->populate(rawdb));
    return std::move(pir_db);
}

Status PIRDatabase::populate(const vector<std::int64_t>& rawdb) {
    if (rawdb.size() != context_->Params()->num_items()) {
        return InvalidArgumentError(
            "Database size " + std::to_string(rawdb.size()) +
            " does not match params value " +
            std::to_string(context_->Params()->num_items()));
    }

    auto evaluator = std::make_unique<seal::Evaluator>(context_->SEALContext());
    db_.resize(rawdb.size());
    for (size_t idx = 0; idx < rawdb.size(); ++idx) {
        try {
            context_->Encoder()->encode(rawdb[idx], db_[idx]);
            if (!context_->Params()->use_ciphertext_multiplication()) {
                evaluator->transform_to_ntt_inplace(
                    db_[idx], context_->SEALContext()->first_parms_id());
            }
        } catch (std::exception& e) {
            return InvalidArgumentError(e.what());
        }
    }
    return absl::OkStatus();
}

Status PIRDatabase::populate(const vector<string>& rawdb) {
    if (rawdb.size() != context_->Params()->num_items()) {
        return InvalidArgumentError(
            "Database size " + std::to_string(rawdb.size()) +
            " does not match params value " +
            std::to_string(context_->Params()->num_items()));
    }

    const auto items_per_pt = context_->Params()->items_per_plaintext();
    db_.resize(context_->Params()->num_pt());
    auto encoder = std::make_unique<StringEncoder>(context_->SEALContext());
    auto evaluator = std::make_unique<seal::Evaluator>(context_->SEALContext());
    if (context_->Params()->bits_per_coeff() > 0) {
        encoder->set_bits_per_coeff(context_->Params()->bits_per_coeff());
    }
    const size_t ITEMS_PER_THREAD = 8000 / items_per_pt * items_per_pt;
    auto rawdb_it = rawdb.begin();
    auto rawdb_end = rawdb.end();
    auto db_size = db_.size();
    const auto use_ctxt_mul = context_->Params()->use_ciphertext_multiplication();
    const auto first_parms_id = context_->SEALContext()->first_parms_id();
    std::vector<std::vector<Status>> ret_status((db_size + ITEMS_PER_THREAD - 1) / ITEMS_PER_THREAD);
    std::vector<std::thread> thrds;
    // for (size_t i = 0; i < (db_size + ITEMS_PER_THREAD - 1) / ITEMS_PER_THREAD; ++i) {
    //     // std::cout << "i = " << i << std::endl;
    //     for (size_t j = 0; j < std::min(ITEMS_PER_THREAD, db_size - i * ITEMS_PER_THREAD); j++) {
    //         // std::cout << "j = " << j << std::endl;
    //         auto end_it = std::min(rawdb_it + items_per_pt, rawdb_end);
    //         ret_status[i].push_back(encoder->encode(rawdb_it, end_it, db_[i * ITEMS_PER_THREAD + j]));
    //         if (!use_ctxt_mul) {
    //             evaluator->transform_to_ntt_inplace(
    //                 db_[i * ITEMS_PER_THREAD + j], first_parms_id);
    //         }
    //         rawdb_it += items_per_pt;
    //     }
    // }

    for (size_t i = 0; i < (db_size + ITEMS_PER_THREAD - 1) / ITEMS_PER_THREAD; ++i) {
        thrds.emplace_back(std::thread(
            [=, &encoder, &evaluator](size_t seg,
                                      std::vector<std::string>::const_iterator curr_begin_it, 
                                      std::vector<Status> &curr_ret_status) {
                for (size_t j = 0; j < std::min(ITEMS_PER_THREAD, db_size - seg * ITEMS_PER_THREAD); j++) {
                    auto end_it = std::min(curr_begin_it + items_per_pt, rawdb_end);
                    curr_ret_status.push_back(encoder->encode(curr_begin_it, end_it, db_[seg * ITEMS_PER_THREAD + j]));
                    if (!use_ctxt_mul) {
                        evaluator->transform_to_ntt_inplace(
                            db_[seg * ITEMS_PER_THREAD + j], first_parms_id);
                    }
                    curr_begin_it += items_per_pt;
                }
            },
            i,
            rawdb_it,
            std::ref(ret_status[i])
            ));
        rawdb_it += items_per_pt * ITEMS_PER_THREAD;
    }
    for (auto &t: thrds) {
        t.join();
    }
    for (const auto& status_vec: ret_status) { 
        for (const auto& status: status_vec) 
            if (!status.ok()) 
                return status; 
    }
    return absl::OkStatus();
}

/**
 * Helper class to make the recursive multiplication operation on the
 * multi-dimensional representation of the database easier. Encapsulates all of
 * the variables needed to do the multiplication, and keeps track of the
 * database iterator to separate it from the database itself.
 */
class DatabaseMultiplier {
   public:
    /**
     * Create a multiplier for the given scenario.
     * @param[in] database Database against which to multiply.
     * @param[in] selection_vector multi-dimensional selection vector
     * @param[in] evaluator Evaluator to use for homomorphic operations.
     * @param[in] relin_keys If not nullptr, relinearization will be done after
     *    every homomorphic multiplication.
     * @param[in] decryptor If not nullptr, outputs to cout the noise budget
     *    remaining after every homomorphic operation.
     */
    DatabaseMultiplier(const vector<Plaintext>& database,
                       vector<Ciphertext>& selection_vector,
                       shared_ptr<Evaluator> evaluator,
                       unique_ptr<CiphertextReencoder> ct_reencoder,
                       std::shared_ptr<seal::SEALContext> seal_context,
                       const seal::RelinKeys* const relin_keys,
                       seal::Decryptor* const decryptor)
        : database_(database),
          selection_vector_(selection_vector),
          evaluator_(evaluator),
          ct_reencoder_(std::move(ct_reencoder)),
          seal_context_(seal_context),
          exp_ratio_(
              ct_reencoder_ == nullptr ? 1 : ct_reencoder_->ExpansionRatio()),
          relin_keys_(relin_keys),
          decryptor_(decryptor) {}

    /**
     * Do the multiplication using the given dimension sizes.
     */
    vector<Ciphertext> multiply(const RepeatedField<uint32_t>& dimensions) {
        database_it_ = database_.begin();
        return multiply(dimensions, selection_vector_.begin(), 0);
    }

   private:
    /**
     * Recursive function to do the dot product of each dimension with the db.
     * Calls itself to move down dimensions until you get to the bottom
     * dimension. Bottom dimension just does a dot product with the DB, and
     * returns the result. Upper levels then take those results, and dot product
     * again with the selection vector, until you get back to the top. NB:
     * Database iterator is kept at the class level so that we move through the
     * database one element at a time.
     *
     * @param[in] dimensions List of remaining demainsion sizes.
     * @param[in] selection_vector_it Iterator into the start of the selection
     *  vector for the current depth.
     * @param[in] depth Current depth.
     */
    vector<Ciphertext> multiply(
        const RepeatedField<uint32_t>& dimensions,
        vector<Ciphertext>::iterator selection_vector_it, size_t depth) {
        const size_t this_dimension = dimensions[0];
        auto remaining_dimensions =
            RepeatedField<uint32_t>(dimensions.begin() + 1, dimensions.end());
        vector<Ciphertext> result;
        auto mul_count = std::min(this_dimension, (size_t)(database_.end() - database_it_));
        if (remaining_dimensions.empty()) {
            // // base case: have to multiply against DB
            // std::vector<std::thread> thrds;
            // std::vector<Ciphertext> partial_results((mul_count + 299) / 300);
            // for (size_t i = 0; i < mul_count; i += 300) {
            //     auto curr_mul_count = std::min(mul_count - i, (size_t)300);
            //     auto &partial_result = partial_results.at(i / 300);
            //     for (size_t j = 0; j < curr_mul_count; j++) {
            //         Ciphertext temp_ct;
            //         if (ct_reencoder_ != nullptr &&
            //             !selection_vector_it[j].is_ntt_form()) {
            //             evaluator_->transform_to_ntt_inplace(
            //                 selection_vector_it[j]);
            //         }
            //         evaluator_->multiply_plain(selection_vector_it[j],
            //                                     database_it_[j], temp_ct);
            //         if (j == 0) {
            //             partial_result = temp_ct;
            //         } else {
            //             evaluator_->add_inplace(partial_result, temp_ct);
            //         }
            //     }
            //     database_it_ += curr_mul_count;
            //     selection_vector_it += curr_mul_count;
            // }

            // base case: have to multiply against DB
            std::vector<std::thread> thrds;
            std::vector<Ciphertext> partial_results((mul_count + 299) / 300);
            for (size_t i = 0; i < mul_count; i += 300) {
                auto curr_mul_count = std::min(mul_count - i, (size_t)300);
                auto &partial_result = partial_results.at(i / 300);
                thrds.emplace_back(std::thread(
                    [=, &partial_result](std::vector<Plaintext>::const_iterator curr_db_it, std::vector<Ciphertext>::iterator curr_select_vec_it) {
                        for (size_t j = 0; j < curr_mul_count; j++) {
                            Ciphertext temp_ct;
                            if (ct_reencoder_ != nullptr &&
                                !curr_select_vec_it[j].is_ntt_form()) {
                                evaluator_->transform_to_ntt_inplace(
                                    curr_select_vec_it[j]);
                            }
                            evaluator_->multiply_plain(curr_select_vec_it[j],
                                                       curr_db_it[j], temp_ct);
                            if (j == 0) {
                                partial_result = temp_ct;
                            } else {
                                evaluator_->add_inplace(partial_result, temp_ct);
                            }
                        }
                    },
                    database_it_,
                    selection_vector_it
                ));
                database_it_ += curr_mul_count;
                selection_vector_it += curr_mul_count;
            }
            for (auto &t: thrds) {
                t.join();
            }

            result.emplace_back(std::move(partial_results[0]));
            for (size_t i = 1; i < partial_results.size(); i++) {
                evaluator_->add_inplace(result[0], partial_results[i]);
            }

            if (ct_reencoder_ != nullptr) {
                evaluator_->transform_from_ntt_inplace(result[0]);
            }
            print_noise(depth, "final", result[0]);
            return result;
        }

        bool first_pass = true;
        for (size_t i = 0; i < this_dimension; ++i) {
            // make sure we don't go past end of DB
            if (database_it_ == database_.end()) break;
            vector<Ciphertext> temp_ct;

            // recurse
            assert(!remaining_dimensions.empty());
            auto lower_result =
                multiply(remaining_dimensions,
                            selection_vector_it + this_dimension, depth + 1);
            print_noise(depth, "recurse", lower_result[0], i);

            // multiply lower results
            if (ct_reencoder_ == nullptr) {
                temp_ct.resize(1);
                evaluator_->multiply(lower_result[0],
                                        *(selection_vector_it + i),
                                        temp_ct[0]);
                print_noise(depth, "mult", temp_ct[0], i);

                if (relin_keys_ != nullptr) {
                    evaluator_->relinearize_inplace(temp_ct[0],
                                                    *relin_keys_);
                    print_noise(depth, "relin", temp_ct[0], i);
                }

            } else {
                // TODO: check that all CT are size 2
                temp_ct.resize(lower_result.size() * exp_ratio_ * 2);
                auto temp_ct_it = temp_ct.begin();
                for (const auto& ct : lower_result) {
                    auto pt_decomp = ct_reencoder_->Encode(ct);
                    size_t k = 0;
                    for (auto pt : pt_decomp) {
                        if (!(selection_vector_it + i)->is_ntt_form()) {
                            evaluator_->transform_to_ntt_inplace(
                                *(selection_vector_it + i));
                        }
                        if (!pt.is_ntt_form()) {
                            evaluator_->transform_to_ntt_inplace(
                                pt, seal_context_->first_parms_id());
                        }
                        evaluator_->multiply_plain(
                            *(selection_vector_it + i), pt, *temp_ct_it);
                        print_noise(depth, "mult", *temp_ct_it, k++);
                        ++temp_ct_it;
                    }
                }
            }

            // add up
            if (first_pass) {
                result = temp_ct;
                first_pass = false;
                print_noise(depth, "first_pass", result[0], i);
            } else {
                for (size_t j = 0; j < result.size(); ++j) {
                    evaluator_->add_inplace(result[j], temp_ct[j]);
                    print_noise(depth, "result", result[j], i);
                }
            }
        }

        for (auto& ct : result) {
            if (ct.is_ntt_form()) {
                evaluator_->transform_from_ntt_inplace(ct);
            }
        }

        print_noise(depth, "final", result[0]);
        return result;
    }

    void print_noise(size_t depth, const string& desc, const Ciphertext& ct,
                     std::optional<size_t> i_opt = {}) {
        if (decryptor_ != nullptr) {
            std::cout << string(depth, ' ');
            if (i_opt) {
                std::cout << "i = " << (*i_opt) << " ";
            }
            std::cout << desc << " noise budget "
                      << decryptor_->invariant_noise_budget(ct) << std::endl;
        }
    }

    const vector<Plaintext>& database_;
    vector<Ciphertext>& selection_vector_;
    shared_ptr<Evaluator> evaluator_;
    unique_ptr<CiphertextReencoder> ct_reencoder_;
    std::shared_ptr<seal::SEALContext> seal_context_;
    const size_t exp_ratio_;

    // If not null, relinearization keys are applied after each HE op
    const seal::RelinKeys* const relin_keys_;

    // If not null, used to get invariant noise budget after each HE op
    seal::Decryptor* const decryptor_;

    // Current location as we move through the database.
    // Needs to be kept here, as lower levels of recursion move forward.
    vector<Plaintext>::const_iterator database_it_;
};

StatusOr<vector<Ciphertext>> PIRDatabase::multiply(
    vector<Ciphertext>& selection_vector,
    const seal::RelinKeys* const relin_keys,
    seal::Decryptor* const decryptor) const {
    auto& dimensions = context_->Params()->dimensions();
    const size_t dim_sum = context_->DimensionsSum();

    if (selection_vector.size() != dim_sum) {
        return InvalidArgumentError(
            "Selection vector size does not match dimensions");
    }

    unique_ptr<CiphertextReencoder> ct_reencoder = nullptr;
    if (!context_->Params()->use_ciphertext_multiplication()) {
        ASSIGN_OR_RETURN(ct_reencoder,
                         CiphertextReencoder::Create(context_->SEALContext()));
    }

    try {
        DatabaseMultiplier dbm(db_, selection_vector, context_->Evaluator(),
                               std::move(ct_reencoder), context_->SEALContext(),
                               relin_keys, decryptor);
        return dbm.multiply(dimensions);
    } catch (std::exception& e) {
        return InternalError(e.what());
    }
}

vector<uint32_t> PIRDatabase::calculate_indices(uint32_t index) {
    uint32_t pt_index = index / context_->Params()->items_per_plaintext();
    vector<uint32_t> results(context_->Params()->dimensions_size(), 0);
    for (int i = results.size() - 1; i >= 0; --i) {
        results[i] = pt_index % context_->Params()->dimensions(i);
        pt_index = pt_index / context_->Params()->dimensions(i);
    }
    return results;
}

size_t PIRDatabase::calculate_item_offset(uint32_t index) {
    uint32_t pt_index = index / context_->Params()->items_per_plaintext();
    return (index - (pt_index * context_->Params()->items_per_plaintext())) *
           context_->Params()->bytes_per_item();
}

vector<uint32_t> PIRDatabase::calculate_dimensions(uint32_t db_size,
                                                   uint32_t num_dimensions) {
    vector<uint32_t> results;
    for (int i = num_dimensions; i > 0; --i) {
        results.push_back(std::ceil(std::pow(db_size, 1.0 / i)));
        db_size = std::ceil(static_cast<double>(db_size) / results.back());
    }
    return results;
}

}  // namespace pir
