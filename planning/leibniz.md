# Leibniz / dy-dx Project — Inventory

Working document. Flat list first, sequencing decided later.
Categories: 
SPINE (must be walked in order, no hand-waving) 
LOAD-BEARING (not on spine, but later spine material won't land without it) 
DIGRESSION (optional, needs its own hook, reader can skip and rejoin with nothing lost).

Thematic tags:
FORMALISES: puts in formal mathematical terms, something that has been introduced or discussed. Given the "no hand waving" rules, ideas are likely to be introduced in a way that avoids unsoundness and that is quite rigorous, but it may lack formal language. But without FORMALISES points made may end up just being assertions, this delivers tightness. This can occur anywhere. It is most likely to occur as a DIGRESSION but it may be part of the SPINE, where a deliverable of the story or point of the story is to explain what something "really is" mathematically.

PORTAL: doesn't need to be rigorous at all, often is closer to sociology/history than mathematics ("some people write it this way because...", "traditional to use k here because..."). Cheap to include, low-stakes, but has real independent payoff: it's what lets someone go and open a real paper. Can sometimes be a part of or the same as FORMALISES but often will not be.

Status tags: idea | rough draft | drafted | needs rewrite | done

---

## SPINE
### Episode 1 - differentials
Headline: "dy/dx is a fraction".
- [ ] Hook: it looks like a fraction and it acts like a fraction. Examples (inversion, chain rule, solving differential equation via differentials to integrate). But we are often told that it isn't. Or that it isn't *rigorously*. That using it is a bit dubious in some unspecified way.
- [ ] Promise: we'll explain what dy and dx "really" are and how you can divide them to get dy/dx.
- [ ] What does `y = x²` actually mean? (is it just an expression in which we insert y's and x's? No, usually we have a picture in mind (picture of it). For this curve we can read off values of y for x. Plotting a graph. Idea of a function. But we can perfectly well write y^2 + x^2 =1 and we get a circle. We can't just read values off. Idea of a constraint. function  — the seed question also sets up the later function/differential split) — status: idea
- [ ] Variables could just be these funny letters to which we assign values in the computer science sense. In the Leibniz sense (not "a number that varies". Start with a curve and see that any point on that curve has numbers associated with it — this is the concept the whole story hangs off. A variable could be defined on the whole plane or just in terms of that single curve. Like s. Still acts as a description of a curve or curves) — status: idea
- [ ] Curves as polygons with infinitesimally small sides. Equivalently as an infinite series of points (the vertices of the polygon). Explanation of why we are not going to use infinitesimals (essential to simplify) but a pointer to an infinitesimal story. Use parametrized curves. Intuition for similarity. Delivery elsewhere. Variables are then functions of t the parameter.
- [ ] If they are functions how do we add/multiply them. Explain pulling back. We can now do any maths we want here. Remark: Very useful idea, used all over.
- [ ] delta of a variable gives a difference. d of a variable gives another function of t. Using limits - refer to load bearing section on this. Now dx is a function also. A differential
- [ ] dy/dx as a genuine ratio of differentials. By way of pulling back we can do maths with dy's and dx's so of course we can take the ratio (the hook / core claim of the whole piece) — status: rough draft
- [ ] The parameter was arbitrary. Observe what happens to it in dy/dx. It vanishes. For equations dy =2dx it doesn't vanish but it doesn't matter (same dependency on both sides). (formal proof here or in a note?).
- [ ] How to make sure that what we end up with has this magic property of not being dependent. Orders of infinity - a dimension system (aside this is a Krull Valuation, how interesting is that and do we need to say much about it here?)
- [ ] Second order. use d on dy/dx with the quotient rule. What you get isn't d^2y/dx^2. What went wrong?
- [ ] If only dx^2 was zero. It would be if it was a constant variable. No reason to assume dx is the same throughout (picture). Leibniz also bothered. "Progression of variables"
- [ ] Link this to parameter dependence. Can we do a parameter change for d^2y/dx^2 and discover the correction factor that way.  this is what makes higher order meaningful later) — status: idea
- [ ] If we pick the progression of variables / parameter we can simplify expressions, but only to a certain extent.
- [ ] That's why the "chain rule" is more complicated/doesn't work. You can't have two progressions being constant.
- [ ] Why dy/dx can be treated just like a fraction (it is really a fraction from a very defensible point of view) but answers may be less useful if you don't respect the orders of infinity system — status: idea

## LOAD-BEARING

- [ ] Delta-method (school) definition, recapped precisely — only introduced
      right before it's needed to contrast against the epsilon/Leibniz approach,
      not at the very start — status: idea
- [ ] Epsilon-style definition, f(x+ε) − f(x) — same deferral logic — status: idea
- [ ] What a differential *is* when there's no obvious function lying around to
      differentiate (this seems to be the crux load-bearing idea — it's the
      thing that makes "variables" necessary rather than decorative) — status: idea
- [ ] Notation: is dy/dx a function of x? How do you indicate applying it —
      dy/dx(1)? dy/dx|_{x=1}? — and why plain dy/dx is fine most of the time
      *because of* the approach being taken (this resolves the friend's
      "what type is y" objection) — status: idea

## DIGRESSION
- [ ] Formalisation: Co-germs, co-germ differentials
- [ ] Why 1-forms don't work (no access to higher derivatives), why forms really don't work (even if integrable, you can make a foliation and then convert all forms to scalars via Hodge duality and then take a d of them, but you need parametrisation for all that - is this true? Can it be fleshed out. Note: it was my original model but it may illustrate why you can't get away from the parameter).
- [ ] Curvature via second order — hook: "what does the *next* term buy you,
      geometrically?". It can illustrate how pleasant it is to play with d. You can just apply d to tan theta=dy/dx and the rest follows very easily algebraically — status: idea
- [ ] A bonus point is that things like ds =dx^2 +dy^2 are natural and in SR you can just do tau = (similar) and derive d tau/dt = gamma (or 1/gamma) directly again just using algebra. All type checks nicely.

- [ ] Hardy on growth-rate infinitesimals (early 20th c. calculus texts defining infinitesimals as classes of function by growth rate) — hook: how do we have infinitely short sides without having infinitesimal numbers? "safe," makes minimal ontological demands  — status: idea
- [ ] Robinson / nonstandard analysis — hook: still needed — status: idea
- [ ] Nilpotents, and Grothendieck's desire for a "double point" — hook: a genuinely motivating picture for why algebraic geometry wants nilpotents at all — status: idea
- [ ] Once we explore "things like ds" we realise we can do inexact differentials, but we didn't know what they were.
- [ ] Explain exactness simply via thing like potentials and voltage. This leads to small-delta notation
- [ ] Small-delta notation for inexact differentials → brief discussion of exactness — hook: still needed — status: idea
- [ ] Bigger digression: thermodynamics. dU=delta q + delta W (how you can guess an exact differential exists if you now something about an inexact one - in this case W along adiabatic paths)
- [ ] Exactness → homology (resist going further than a taste) — hook: still
      needed, and a place to explicitly signal "you can stop here" — status: idea

## OPEN QUESTIONS (mine, not the reader's)

- Does the dy/dx(1) notation discussion belong in LOAD-BEARING or is it
  actually a DIGRESSION that just feels urgent because it bugged a real person?
- Where exactly does "variables in the Leibniz sense" first get said —
  before or after the y=x² seed question?
- Is curvature-via-second-order really a digression, or is it secretly
  load-bearing for something later (e.g. if parametrised curves / speed in
  the tangent bundle ever enters, à la Burke)?
- How much of the infinitesimals material (Hardy / Robinson / nilpotents) is
  one digression with three branches, vs three separate digressions? I wonder if there isn't an "infinitesimals" story - a sort of episode three, which can be motivated directly by the question of how is a curve a polygon. OR an alternative is that we can insert the function growth rate definition very early because it so natural; was widely used; and as a digression or further digression it is the basis of lots of O(x) type notation.

## PARKING LOT (unsorted — anything new goes here first, sort later)

-

