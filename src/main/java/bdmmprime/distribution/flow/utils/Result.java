/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.distribution.flow.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Stream;

/**
 * This is a wrapper class which holds the result of a computation or a potential error
 * which was thrown during that computation.
 * This is useful when running things in parallel, as it allows to throw potential errors
 * on the main sequential thread once all threads have been completed. If the threads themselves
 * were to throw the errors directly, the other threads would finish in the background leading
 * to potential deadlocks.
 */
public class Result<T> {

    private final T result;
    private final RuntimeException error;

    public Result(T result, RuntimeException error) {
        this.result = result;
        this.error = error;
    }

    /**
     * Executes the supplier and returns a Result object containing the result or the error thrown.
     * Any error is caught and stored in the returned object instead of thrown.
     */
    public static <T> Result<T> of(Supplier<T> supplier) {
        try {
            return Result.success(supplier.get());
        } catch (RuntimeException error) {
            return Result.failure(error);
        }
    }

    public static <T> Result<T> success(T result) {
        return new Result<>(result, null);
    }

    public static <T> Result<T> failure(RuntimeException error) {
        return new Result<>(null, error);
    }

    /** Consumes the given stream of result objects and throws the first error encountered. */
    public static <T> void throwIfFailure(Stream<Result<T>> results) {
        Result.getAllOrThrow(results);
    }

    /**
     * Consumes the given stream of result objects and returns the unwrapped results. The stream is
     * fully evaluated before the first error encountered is thrown, so no computation of a parallel
     * stream keeps running in the background.
     */
    public static <T> List<T> getAllOrThrow(Stream<Result<T>> results) {
        List<Result<T>> resultList = results.toList();

        List<T> unwrappedResults = new ArrayList<>(resultList.size());
        for (Result<T> result : resultList) {
            unwrappedResults.add(result.getOrThrow());
        }
        return unwrappedResults;
    }

    /** Unwraps the stored result or throws the stored error. */
    public T getOrThrow() {
        if (this.error != null) {
            throw this.error;
        }
        return this.result;
    }
}
