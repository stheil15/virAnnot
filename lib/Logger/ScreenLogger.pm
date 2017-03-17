#!/share/apps/bin/perl
# $Id: Logger.pm,v 1.3 2011/07/01 18:07:31 nguilhot Exp $
package Logger::ScreenLogger;

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
use Exporter 'import';

our $logger;
our @EXPORT = qw($logger);

$logger = Log::Log4perl->get_logger('');

# Define stdout Appender, by default messages will only be logged to stdout
# In order to log messages to a file, you need to call the initFileLoggers
# method with the name of folder where to create the log files
my $stdout_layout = Log::Log4perl::Layout::PatternLayout->new("%5p - %m%n");
my $stdout_appender =  Log::Log4perl::Appender->new(
        "Log::Log4perl::Appender::Screen",
        name      => 'screenlog',
        stderr    => 1);
$stdout_appender->layout($stdout_layout);
$stdout_appender->threshold($INFO);
$logger->add_appender($stdout_appender);
$logger->level($INFO);



sub changeMode {
    my ($verbosity) = @_;
    my $debug_layout = Log::Log4perl::Layout::PatternLayout->new("%d %5p> %F{1}:%L %M - %m%n");
    my $debug_appender = Log::Log4perl::Appender->new(
        "Log::Log4perl::Appender::Screen",
        name      => 'debuglog',
        stderr    => 0);
    $debug_appender->layout($debug_layout);
    $logger->add_appender($debug_appender);
    $logger->remove_appender('screenlog');

    if($verbosity == 3){
        $stdout_appender->threshold($DEBUG);
        $logger->level($DEBUG);
    }
    if($verbosity == 4){

        $stdout_appender->threshold($TRACE);
        $logger->level($TRACE);
    }
}

1;
