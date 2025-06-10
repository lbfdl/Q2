function test_fLJsum()
%TEST_FLJSUM Simple verification of forces computed by fLJsum
%   Compares force results from fLJsum with the analytical
%   Lennard-Jones expression for two-particle systems.

    a = 10;      % simulation box length
    rc = a/2;    % cutoff large enough for these pairs
    n = 2;

    %% First pair: r = 1.0 (repulsive)
    r = 1.0;
    x = [0; r/a];
    y = [0; 0];
    z = [0; 0];
    [~, ~, fx, ~, ~] = fLJsum(a, n, rc, x, y, z);
    Fexpected = -24*(2/r^13 - 1/r^7);
    assert(abs(fx(1) - Fexpected) < 1e-6, 'Force mismatch for r=1.0');
    assert(abs(fx(2) + Fexpected) < 1e-6, 'Newton''s third law violated for r=1.0');

    %% Second pair: r = 1.5 (attractive)
    r = 1.5;
    x = [0; r/a];
    [~, ~, fx, ~, ~] = fLJsum(a, n, rc, x, y, z);
    Fexpected = -24*(2/r^13 - 1/r^7);
    assert(abs(fx(1) - Fexpected) < 1e-6, 'Force mismatch for r=1.5');
    assert(abs(fx(2) + Fexpected) < 1e-6, 'Newton''s third law violated for r=1.5');

    disp('All tests passed.');
end
