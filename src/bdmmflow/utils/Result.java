package bdmmflow.utils;

import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;
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
    T result;
    RuntimeException error;

    public Result(T result, RuntimeException error) {
        this.result = result;
        this.error = error;
    }

    /**
     * Executes the supplier and returns a Result error containing the result or error thrown.
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
        List<Result<T>> resultList = results.toList();
        for (Result<T> result : resultList) {
            result.getOrThrow();
        }
    }

    /** Unwraps the stored result or throws the stored error. */
    public T getOrThrow() {
        if (this.error != null) {
            throw this.error;
        }
        return this.result;
    }
}
