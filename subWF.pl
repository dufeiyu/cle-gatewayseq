#! /usr/bin/perl

#Copyright (C) 2021 Feiyu Du <fdu@wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;

umask 002;

use Cwd qw(abs_path);
use JSON qw(from_json to_json);
use IO::File;
use File::Spec;
use File::Basename;

die "Provide GatewaySeq output directory" unless @ARGV == 1;
my $dir = $ARGV[0];
#my $dir = '/storage1/fs1/duncavagee/Active/SEQ/MyeloSeqHD/batchdir/CI-733_rerun';
die "$dir is not a valid directory" unless -d $dir;
$dir = abs_path($dir);

my $git_dir = '/storage1/fs1/duncavagee/Active/SEQ/GatewaySeq/process/git/cle-gatewayseq';
my $conf = File::Spec->join($git_dir, 'application.conf');
my $wdl  = File::Spec->join($git_dir, 'GatewayseqAnalysis.wdl');
my $json_template = File::Spec->join($git_dir, 'GatewayseqAnalysis.json');

my $group  = '/cle/wdl/haloplex';
my $queue  = 'dspencer';
my $docker = 'registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-38';

my $user_group = 'compute-gtac-mgi';

my $ct = 0;
for my $case_dir (glob("$dir/TW*"), glob("$dir/H_*")) {
    my ($case_name) = basename($case_dir) =~ /^(\S+)_[ATCG]{10}\-/;
    die "Fail to get case name for $case_dir" unless $case_name;
    my $dragen_dir = File::Spec->join($case_dir, 'dragen');
    die "$dragen_dir not existing" unless -d $dragen_dir;
    
    my $cram = File::Spec->join($dragen_dir, $case_name.'_tumor.cram');
    my $cram_index = $cram.'.crai';
    my $vcf = File::Spec->join($dragen_dir, $case_name.'.hard-filtered.vcf.gz');
    my $vcf_index = $vcf.'.tbi';
    my $svvcf = File::Spec->join($dragen_dir, $case_name.'.sv.vcf.gz');
    my $svvcf_index = $vcf.'.tbi';

    die "one, some, or all of cram, vcf, svvcf and their index not valid" unless -s $cram and -s $cram_index and -s $vcf and -s $vcf_index and -s $svvcf and -s $svvcf_index;

    my $inputs = from_json(`cat $json_template`);
    $inputs->{'GatewayseqAnalysis.DragenCram'}       = $cram;
    $inputs->{'GatewayseqAnalysis.DragenCramIndex'}  = $cram_index;
    $inputs->{'GatewayseqAnalysis.DragenVcf'}        = $vcf;
    $inputs->{'GatewayseqAnalysis.DragenVcfIndex'}   = $vcf_index;
    $inputs->{'GatewayseqAnalysis.DragenSvVcf'}      = $svvcf;
    $inputs->{'GatewayseqAnalysis.DragenSvVcfIndex'} = $svvcf_index;

    $inputs->{'GatewayseqAnalysis.Name'}             = $case_name;
    $inputs->{'GatewayseqAnalysis.OutputDir'}        = $dir;
    $inputs->{'GatewayseqAnalysis.SubDir'}           = basename($case_dir);

    my $input_json = File::Spec->join($case_dir, 'GatewayseqAnalysis.json');
    my $json_fh = IO::File->new(">$input_json") or die "fail to write to $input_json";

    $json_fh->print(to_json($inputs, {canonical => 1, pretty => 1}));
    $json_fh->close;

    my $out_log = File::Spec->join($case_dir, 'out.log');
    my $err_log = File::Spec->join($case_dir, 'err.log');

    my $cmd = "bsub -g $group -G $user_group -oo $out_log -eo $err_log -q $queue -R \"select[mem>8000] rusage[mem=8000]\" -M 8000000 -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl -i $input_json $wdl";

    system $cmd;
    #print $cmd."\n";
    print $case_name." submitted\n";
    $ct++;
    #last if $ct == 1;
    sleep 10; #DB upload and query
}
print "All $ct done\n";
