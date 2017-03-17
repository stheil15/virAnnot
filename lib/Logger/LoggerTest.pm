#!/usr/bin/perl

package Logger::LoggerTest;

##################################################
## Included modules
##################################################
## Perl modules
use strict;
use warnings;
use diagnostics;
use Log::Log4perl;
use Log::Log4perl::Layout;
use Log::Log4perl::Level;
use Log::Log4perl::Filter::LevelRange;
use Exporter 'import';

our $logger;
our @EXPORT = qw($logger);

$logger = Log::Log4perl->get_logger('');

foreach my $appender (keys Log::Log4perl->appenders()){

	Log::Log4perl->eradicate_appender($appender);
}

############################## CONFIGURATION OF STDOUT #################################

my $stdout_layout = Log::Log4perl::Layout::PatternLayout->new("%5p - %m%n");

my $stdout_filter = Log::Log4perl::Filter::LevelRange->new(
    LevelMin      => 'ALL',
    LevelMax      => 'INFO',
    AcceptOnMatch    => 1
);

############################## CONFIGURATION OF STDERR #################################

my $stderr_layout = Log::Log4perl::Layout::PatternLayout->new("%5p - %m%n");

my $stderr_filter = Log::Log4perl::Filter::LevelRange->new(
    LevelMin      => 'WARN',
    LevelMax      => 'FATAL',
    AcceptOnMatch    => 1
);

############################## CONFIGURATION OF LOGGER #################################

$logger->level('ALL');

Logger::LoggerTest->changeMode(3);
Logger::LoggerTest->setErrorCutoff(2);

Logger::LoggerTest->setInfoOutput();
Logger::LoggerTest->setErrorOutput();

############################### SOME TESTS OF LOGGER ###################################

#Logger::LoggerTest->setInfoOutput('info.log');
#Logger::LoggerTest->setErrorOutput('error.log');

#$logger->trace("Here's the trace (5)");
#$logger->debug("Here's the debug (4)");
#$logger->info("Here's the info (3)");
#$logger->warn("Here's the warn (2)");
#$logger->error("Here's the error (1)");
#$logger->fatal("Here's the fatal (0)");

#Logger::LoggerTest->setInfoOutput();
#Logger::LoggerTest->setErrorOutput();

#$logger->trace("Here's other trace (5)");
#$logger->debug("Here's other debug (4)");
#$logger->info("Here's other info (3)");
#$logger->warn("Here's other warn (2)");
#$logger->error("Here's other error (1)");
#$logger->fatal("Here's other fatal (0)");

=head1 LOGGER METHODS

=head2

=head2 getVerbosityLevels

=head2

=head3 Arguments

=over 4

=item

None.

=back

=head3 Returns

=over 4

=item

A hash where keys are numeric verbosity levels and values the corresponding name of levels in Log4perl

=back

=cut

sub getVerbosityLevels{

    my %verbosityLevels = (

        '0' =>  'FATAL',
        '1' =>  'ERROR',
        '2' =>  'WARN',
        '3' =>  'INFO',
        '4' =>  'DEBUG',
        '5' =>  'TRACE',
    );

    return %verbosityLevels
}

=head2

=head2 getLogFilesMode

=head2

=head3 Arguments

=over 4

=item

None.

=back

=head3 Returns

=over 4

=item

A hash where keys are available mode for opening a log file and values the corresponding mode names in Log4perl

=back

=cut

sub getLogFilesMode{

    my %logFilesMode = (

        'overwrite' =>  'write',
        'write'     =>  'write',
        'preserve'  =>  'append',
        'append'    =>  'append',
        'truncate'  =>  'truncate',
    );

    return %logFilesMode
}

=head2

=head2 changeMode

=head2

=head3 Arguments

=over 4

=item

A numeric verbosity level between 0 and 5 (see getVerbosityLevels for more information) to apply to the logger

=back

=head3 Returns

=over 4

=item

None.

=back

=cut

sub changeMode {

    my ($class, $verbosity) = @_;

    my %verbosityLevels = $class->getVerbosityLevels;

    if ( defined $verbosity ){

            if ( exists $verbosityLevels{$verbosity} ){

                $logger->level( $verbosityLevels{$verbosity} );
            }

            else{

               $logger->warn('Verbosity level must be a number between 0 and 5');
            }
    }
}

=head2

=head2 setErrorCutoff

=head2

=head3 Arguments

=over 4

=item

An optional numeric verbosity level cutoff between 0 and 5 (see getVerbosityLevels for more information) for error messages.

All the message with a verbosity level >= verbosity cutoff will be considered as error messages.

If no verbosity level is provided, none of the messages will be considered as error.

=back

=head3 Returns

=over 4

=item

None.

=back

=cut

sub setErrorCutoff {

    my ($class, $verbosity) = @_;

    my %verbosityLevels = $class->getVerbosityLevels;

    if ( defined $verbosity ){

            if ( exists $verbosityLevels{$verbosity} ){

                $stderr_filter->{'LevelMin'} = $verbosityLevels{$verbosity};

                if ( exists $verbosityLevels{$verbosity+1} ){

                    $stdout_filter->{'LevelMax'} = $verbosityLevels{$verbosity+1};
                }

                else{

                    $stdout_filter->{'LevelMax'} = 'ALL';
                }
            }

            else{

                $logger->warn('Verbosity cutoff level must be a number between 0 and 5');
            }
    }

    else{

        $stderr_filter->{'LevelMin'} = 'OFF';
        $stdout_filter->{'LevelMax'} = 'FATAL';
    }
}

=head2

=head2 setErrorOutput

=head2

=head3 Arguments

=over 4

=item

An output file where to write error messages (optional). If no defined, the error messages will be written on STDERR

=item

An opening mode for the log file (optional) :

- overwrite | write : the new logs will overwritte the previous logs (equivalent to >).

- preserve | append : the new logs will be written at the end of the file without deleting the previous log (equivalent to >>)

- truncate : the log file will be truncated

see getLogFilesMode for more information

=back

=head3 Returns

=over 4

=item

None.

=back

=cut

sub setErrorOutput{

    my ($class, $output, $write_mode) = @_;

    my $stderr_appender;
    my $name = 'errorlog';

    if(Log::Log4perl->appender_by_name($name)){

        $logger->remove_appender($name);
    }

    if(defined $output){

        my %logFilesMode = $class->getLogFilesMode;

        if(! defined $write_mode){

            $write_mode = 'preserve';
        }

        $write_mode = $logFilesMode{$write_mode};

        $stderr_appender =  Log::Log4perl::Appender->new(
            "Log::Log4perl::Appender::File",
            name      => $name,
            filename  => $output,
            mode      => $write_mode,
        );
    }

    else{

        $stderr_appender =  Log::Log4perl::Appender->new(
            "Log::Log4perl::Appender::Screen",
            name      => $name,
            stderr    => 1,
        );
    }

    $stderr_appender->layout($stderr_layout);
    $stderr_appender->filter($stderr_filter);
    $logger->add_appender($stderr_appender);
}

=head2

=head2 setInfoOutput

=head2

=head3 Arguments

=over 4

=item

An output file where to write information and debug messages (optional). If no defined, the information and debug messages will be written on STDOUT

=item

An opening mode for the log file (optional) :

- overwrite | write : the new logs will overwritte the previous logs (equivalent to >).

- preserve | append : the new logs will be written at the end of the file without deleting the previous log (equivalent to >>)

- truncate : the log file will be truncated

see getLogFilesMode for more information

=back

=head3 Returns

=over 4

=item

None.

=back

=cut

sub setInfoOutput{

    my ($class, $output, $write_mode) = @_;

    my $stdout_appender;
    my $name = 'infolog';

    if(Log::Log4perl->appender_by_name($name)){

        $logger->remove_appender($name);
    }

    if(defined $output){

        my %logFilesMode = $class->getLogFilesMode;

        if(! defined $write_mode){

            $write_mode = 'preserve';
        }

        $write_mode = $logFilesMode{$write_mode};

        $stdout_appender =  Log::Log4perl::Appender->new(
            "Log::Log4perl::Appender::File",
            name      => $name,
            filename  => $output,
            mode      => $write_mode,
        );
    }

    else{

        $stdout_appender =  Log::Log4perl::Appender->new(
            "Log::Log4perl::Appender::Screen",
            name      => $name,
            stderr    => 0,
        );
    }

    $stdout_appender->layout($stdout_layout);
    $stdout_appender->filter($stdout_filter);
    $logger->add_appender($stdout_appender);
}

1;
