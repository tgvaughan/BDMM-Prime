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
