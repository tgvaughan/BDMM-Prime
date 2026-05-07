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

package bdmmprime.util;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Logger;

import java.io.IOException;

@Description("A Logger which can be disabled in BEAUti via a checkbox.")
public class OptionalLogger extends Logger {

    public Input<Boolean> enableLoggerInput = new Input<>("enableLogger",
            "If false, logger produces no output. Default is false.", false);

    @Override
    public void initAndValidate() {
        if (enableLoggerInput.get())
            super.initAndValidate();
    }

    @Override
    public void init() throws IOException {
        if (enableLoggerInput.get())
            super.init();
    }

    @Override
    public void log(long sampleNr) {
        if (enableLoggerInput.get())
            super.log(sampleNr);
    }

    @Override
    public void close() {
        if (enableLoggerInput.get())
            super.close();
    }

    @Override
    public boolean isLoggingToStdout() {

        // In the case that the logger is not enabled, report that we're
        // logging to stdout.  This is a hack to avoid problems with
        // missing files and state resumption.

        if (enableLoggerInput.get())
            return super.isLoggingToStdout();
        else
            return true;
    }
}
