#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 2) {
    die "Usage: perl trop_v7.pl <trf_output_file> <sequence_lengths_file>\n";
}

my ($trf_file, $seq_lengths_file) = @ARGV;

# Step 1: Read sequence lengths
print "Importing sequence lengths...\n";
my %seq_lengths;
open my $seq_fh, '<', $seq_lengths_file or die "Cannot open $seq_lengths_file: $!";
while (<$seq_fh>) {
    chomp;
    my @cols = split(/\t/, $_);
    $seq_lengths{$cols[0]} = $cols[1] if @cols == 2;
}
close $seq_fh;
print "Found ", scalar(keys %seq_lengths), " sequences\n";

# Step 2: Read TRF output and parse data
print "Importing and parsing TRF file...\n";
my %trf_data;
my $current_seq = "";
open my $trf_fh, '<', $trf_file or die "Cannot open $trf_file: $!";
while (<$trf_fh>) {
    chomp;
    if (/^Sequence:\s*(\S+)/) {  
        $current_seq = $1;  
    } elsif (/^\d+\s+\d+/) {  
        my @cols = split(/\s+/, $_);
        if (defined $current_seq) {
            push @{$trf_data{$current_seq}}, [$cols[0], $cols[1]];
        }
    }
}
close $trf_fh;
print "Found tandems from ", scalar(keys %trf_data), " sequences\n";

# Step 3: Merge overlapping tandems
print "Merging overlapping tandems...\n";
my %merged_trf;
foreach my $seq (keys %trf_data) {
    my @sorted = sort { $a->[0] <=> $b->[0] } @{$trf_data{$seq}};
    my @merged;
    my ($minpos, $maxpos) = @{$sorted[0]};  

    for (my $i = 1; $i < @sorted; $i++) {
        my ($start, $end) = @{$sorted[$i]};
        if ($start <= $maxpos) {  
            $maxpos = $end if $end > $maxpos;
        } else {  
            push @merged, [$minpos, $maxpos];
            ($minpos, $maxpos) = ($start, $end);
        }
    }
    push @merged, [$minpos, $maxpos];
    $merged_trf{$seq} = \@merged;
}

# Step 4: Export combined file (1_tandem_per_read.txt)
open my $out1, '>', "1_tandem_per_read.txt" or die "Cannot write file\n";
print $out1 "seqname\ttrgroup\ttrlength\tseqlength\ttrprop\tstart_5prime\tstop_5prime\tstart_3prime\tstop_3prime\n";

my %tr_summary;
foreach my $seq (sort keys %merged_trf) {
    my $seq_len = $seq_lengths{$seq} // 0;
    my $trgroup = 0;
    my $total_trlength = 0;

    foreach my $tr (@{$merged_trf{$seq}}) {
        my ($start, $end) = @$tr;
        my $tr_length = $end - $start;
        $total_trlength += $tr_length;
        $trgroup++;

        my $start_3prime = ($seq_len - $end) * (-1) - 1;
        my $stop_3prime  = ($seq_len - $start) * (-1) - 1;
        my $tr_prop = $seq_len ? sprintf("%.2f", ($tr_length * 100) / $seq_len) : 0;

        print $out1 "$seq\t$trgroup\t$tr_length\t$seq_len\t$tr_prop\t$start\t$end\t$start_3prime\t$stop_3prime\n";
    }

    $tr_summary{$seq} = { total => $total_trlength, trgroup => $trgroup };
}
close $out1;
print "File 1 exported\n";

# Step 5: Compute summary statistics (2_tandem_summary.txt)
open my $out2, '>', "2_tandem_summary.txt" or die "Cannot write file\n";
print $out2 "class\tnseq\tnseq_perc\tbp\tbp_perc\n";

my $total_reads = scalar(keys %seq_lengths);
my $total_bp = 0;
my $tandem_bp = 0;
my $tandem_reads = 0;

foreach my $seq (keys %seq_lengths) {
    $total_bp += $seq_lengths{$seq};
    if (exists $tr_summary{$seq}) {
        $tandem_bp += $tr_summary{$seq}{total};
        $tandem_reads++;
    }
}

my $no_tandem_reads = $total_reads - $tandem_reads;
my $no_tandem_bp = $total_bp - $tandem_bp;

print $out2 "No Tandem\t$no_tandem_reads\t", sprintf("%.2f", ($no_tandem_reads * 100) / $total_reads), "\t$no_tandem_bp\t", sprintf("%.2f", ($no_tandem_bp * 100) / $total_bp), "\n";
print $out2 "Tandem\t$tandem_reads\t", sprintf("%.2f", ($tandem_reads * 100) / $total_reads), "\t$tandem_bp\t", sprintf("%.2f", ($tandem_bp * 100) / $total_bp), "\n";
print $out2 "Total\t$total_reads\t100.00\t$total_bp\t100.00\n";

close $out2;
print "File 2 exported\n";
print "Processing completed.\n";

