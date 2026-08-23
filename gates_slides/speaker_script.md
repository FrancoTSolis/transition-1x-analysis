# Speaker script — Gates Foundation visit

*Two slides, about 3 minutes. Plain English. Click or press → to advance.*

---

## Slide 1 — AI–Quantum Co-Design (~85 s)

> Batteries, catalysts, medicines — they all come down to the same question:
> what are the molecules doing? Simulating that accurately is the thing
> standing between us and a lot of new technology.
>
> There are two tools for it, and each has a hole in it.
>
> **[left]** A quantum computer speaks nature's own language. Point it at one
> molecule and it gives a near-exact answer, including the chemistry ordinary
> supercomputers get wrong. But it is one molecule at a time, on machines that
> stay small and expensive.
>
> **[right]** AI is the mirror image. It handles thousands of molecules at
> once, instantly — but it is only as good as the data it learned from, and
> most of that data today comes from approximate methods.
>
> **[the loop]** What we are building is not one or the other; it is the loop
> between them. Across the top: the quantum computer's exact answers become
> training data the AI can trust. Around the bottom: the AI designs and
> configures the next quantum run, so the machine is spent only where it
> matters.
>
> That is what co-design means here — the two are developed against each
> other, not in sequence.
>
> **[bottom bar]** And it compounds: better training data, machines that are
> actually usable, bigger molecules within reach, and every exact answer
> producing a better next model.

---

## Slide 2 — Co-design in practice: chemistry on a quantum computer (~90 s)

> Here is one turn of that loop, already built and measured — and it shows how
> a quantum computer actually does chemistry.
>
> **[left]** You start with a molecule.
>
> **[middle, the circuit]** A quantum computer doesn't "run chemistry"
> directly. You have to build a circuit for that specific molecule. The one we
> use is called the LUCJ ansatz, and at a high level it alternates two kinds of
> operation: the blue blocks rotate the electron orbitals, the orange blocks
> let the electrons interact — repeated in layers.
>
> The catch is that every one of those boxes is filled with numbers, hundreds
> of them, and they are different for every molecule. Finding them is a
> separate optimization done from scratch each time — hours of computing before
> any chemistry happens.
>
> **[the gold box]** That is the step our AI replaces. It reads the molecule
> and writes all those numbers in a single pass.
>
> And this is where the real saving is. The numbers decide where the circuit
> starts. Start far from the answer and the machine needs many rounds of
> running, measuring and adjusting, and every one of those rounds costs
> quantum time. Start close, and it settles in far fewer runs. So the AI is
> not just saving classical setup time, it is buying back time on the quantum
> computer itself.
>
> **[right]** Then the circuit runs on quantum hardware and gives back the
> molecular energy — the number chemistry actually needs.
>
> **[the dashed return arrow]** And here is where the loop from the first
> slide closes. That energy is also a score: it tells us how good the AI's
> circuit actually was. Feed it back as a reward and the model tunes itself
> against the quantum computer's own answers. That is the step we are building
> now.
>
> **[right panel]** As for how well it works today: each dot is one molecule
> the model had never seen — the exact value across, the AI's prediction up.
> Perfect prediction means sitting on the line, and they do: R-squared 0.96.

**Closing line:** "AI makes quantum computers practical; quantum computers
make AI trustworthy. We build them as one system."

---

## Short answers if asked

**Is the whole loop running?** The AI-designs-the-quantum-run half is working
at the scale on slide 2. Feeding quantum results back as the training signal
is what we are building next.

**Does 91% mean 91% accurate energies?** No — it is how closely the predicted
configuration matches the exact one. Energy-level validation is the next
experiment, and we are careful not to claim it yet.

**Why not just build a bigger quantum computer?** Even a big one is wasted if
every run needs hours of setup. And the same idea lets a small quantum machine
handle just the hard fragment of a much larger molecule.

**What would help most?** Running this on real quantum hardware rather than
simulators, and releasing the datasets openly.
