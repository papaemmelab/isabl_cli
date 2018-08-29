class AbstractBatchSystem():

    def submit_analyses(self, app, command_tuples):  # pylint: disable=W9008
        """
        Submit analyses.

        Arguments:
            app (cli.AbstractApplication): app instance.
            command_tuples (list): list of (analysis, command) tuples.

        Returns:
            list: of analysis, status tuples.
        """
        raise NotImplementedError
