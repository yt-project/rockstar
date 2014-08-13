#!/usr/bin/perl -w
print "Enter number of pipes to create (usually, the same as the number\n of snapshots): ";
chomp(my $num_pipes = <>);
print "Enter path to directory where pipes should be created [default: .]: ";
chomp(my $dir = <>);
$dir = "." unless (length $dir);

for (1..$num_pipes) {
    my $pfn = sprintf("%s/pipe.%d", $dir, $_-1);
    if (system("mkfifo", $pfn) != 0) {
	die "Error creating named pipe $pfn!\n";
    } else {
	print "Created $pfn.\n";
    }
}
