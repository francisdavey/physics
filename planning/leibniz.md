# Project: Leibniz / dy-dx Project

## Document structure
Working document. Flat list first, sequencing decided later. the list has a slightly hierarchical structure, using ### section tags. These are blobs of thought that probably go together within any category. Order of sections is not decided and membership of a section is a soft grouping idea.

Categories:  
SPINE (must be walked in order, no hand-waving). One test for this is *if you removed this item, does the piece still deliver on its central promise, or does the promise go unfulfilled?* (per Claude)  
LOAD-BEARING (not on spine, but later spine material won't land without it)  
DIGRESSION (optional, needs its own hook, reader can skip and rejoin with nothing lost).  
PARKING LOT (very optional, for things I have discovered or thought about, that don't have any obvious place. They may be seeds for new articles or just reference for myself).

Thematic tags:
FORMALISES: puts in formal mathematical terms, something that has been introduced or discussed. Given the "no hand waving" rules, ideas are likely to be introduced in a way that avoids unsoundness and that is quite rigorous, but it may lack formal language. But without FORMALISES points made may end up just being assertions, this delivers tightness. This can occur anywhere. It is most likely to occur as a DIGRESSION but it may be part of the SPINE, where a deliverable of the story or point of the story is to explain what something "really is" mathematically.    

PORTAL: doesn't need to be rigorous at all, often is closer to sociology/history than mathematics ("some people write it this way because...", "traditional to use k here because..."). Cheap to include, low-stakes, but has real independent payoff: it's what lets someone go and open a real paper. Can sometimes be a part of or the same as FORMALISES but often will not be.    

Sketch: is a tag for some text that means I am expanding out an idea. This may be used in the final write-up, but helps feel out what is in my head.

REMARK: tags tell are things said to myself. Sometimes these will include remarks on formalisation

Caveat: is a place where what I am saying is not quite true. Those need to be flagged and paid off. Eg, if I say "functions have only one value" that is true - in a aense - for a definition of "function" but "multi-valued function" is often used. That can be formalised.

---
## SPINE
### Episode  - differentials
Headline: "dy/dx is a fraction".  
-  Hook: it looks like a fraction and it acts like a fraction. Examples (inversion, chain rule, solving differential equation via differentials to integrate). But we are often told that it isn't. Or that it isn't *rigorously* a fraction. That using it is a bit dubious in some unspecified way.
-  Promise: we'll explain what dy and dx "really" are and how you can divide them to get dy/dx.
-  Aims of piece, to give a thorough account of how differentials may be used "in the spirit of Leibniz". Two strategic choices seem sensible: (1) is to give the account in terms of limits rather than infinitesimals; postponing infinitesimals to a follow-on series or post. It may turn out to be sensible to touch on infinitesimals a little, that remains to be seen. (2) Integration above 1D seems to introduce considerable further questions, which may require a separate or careful approach.

#### tell them what dx etc means.
- What do equations mean. Eg what does `y = x²` actually mean?
-- Sketch (1) Is it just an expression in which we insert y's and x's? No, usually we have a picture in mind. PICTURE (graph of equation) For this curve we can read off values of y for x. This is essentially what a "function" is, when our teachers said "y is a function of x" they had this in mind. "Plotting a graph". 
-- (2) But we can perfectly well write y^2 + x^2 =1 and we get a circle. We can't just read values off. No function here [Caveat: partial functions are an idea as are "multi-valued" functions]
-- (3) Important to us because Leibniz wasn't thinking about functions. That came later.
-- (4) Better understood as a constraint or set of conditions a curve must follow.
- "Variables" (segue: before we can  talk about what dx and dy might be, we have to answer what x and y are. They obviously aren't just stand-ins for numbers because we can't meaningfully say one number "depends on" another, let alone is a "function of" another,  - answer: they are "variables"). Introduce Leibniz's idea of a variable as something associated with a point on a curve.
-- sketch: (1) Variables could just be these funny letters to which we assign values in the computer science sense.
  (2) Leibniz thought of them as things associated with points on a curve. They varied with the curve.
-- (3) possibly: Draw a curve with variables associated with it. Perhaps contrast x and y (functions on the plane) with s (only defined on a curve).
- Define how to perform arithmetic on variables. PORTAL:  idea of "pulling back", popular and  much used.
- Curves as polygons with infinitesimally small sides. [remark: could we comment on the proof of the area of a circle being pi r squred, this was how my elementary school teacher proved it) Equivalently as an infinite series of points (the vertices of the polygon).
-  Explanation of why we are not going to use infinitesimals (essential to simplify) but a pointer to an infinitesimal story.
-  Use parametrized curves instead of infinitesimals. Variables are then functions of t the parameter. It is much easier to present the theory as x(t) dx(t) than in any other way, so we need this marker somewhere. Parametrization is really much the same as having really small sides we are just saying that the reals do the job.
- Delta of a variable gives a difference. Since we are thinking of curves as polygons, delta is well defined.
- Transition to "dx" in the non-infinitesimal world. It computes x(t + epsilon) -x over epsilon in the limit. This is another function of t. We are using limits - refer to load bearing section on this. Now dx is a function also.
- Can repeat for df for any "variable".
- Name: a differential
- We can do arithmetic (any real maths we like) on differentials (since they are variables) via "pulling back" again.
- Pulling it together: dy/dx as a genuine ratio of differentials. By way of pulling back we can do maths with dy's and dx's so of course we can take the ratio (the hook / core claim of the whole piece)

####  parameter dependence and order of infinity - how to make sure your expressions are meaningful
- [ ] The parameter was arbitrary. Observe what happens to it in dy/dx. It vanishes. For equations dy =2dx it doesn't vanish but it doesn't matter (same dependency on both sides). (formal proof here or in a note?).
- [ ] How to make sure that what we end up with has this magic property of not being dependent. Orders of infinity - a dimension system (PORTAL to ? )
- [ ] Caveats with order of infinity - constrains the kind of expression we can usefully work with, eg log(dx).

#### Second order
- Second order. use d on dy/dx with the quotient rule. What you get isn't d^2y/dx^2. What went wrong?
- If only dx^2 was zero. It would be if it was a constant variable. No reason to assume dx is the same throughout (picture). Leibniz also bothered. "Progression of variables"
- Link this to parameter dependence. Can we do a parameter change for d^2y/dx^2 and discover the correction factor that way.  this is what makes higher order meaningful later) 
-  If we pick the progression of variables / parameter we can simplify expressions, but only to a certain extent.
-  That's why the "chain rule" is more complicated/doesn't work. You can't have two progressions being constant.
- Conclusion: Why dy/dx can be treated just like a fraction (it is really a fraction from a very defensible point of view) but answers may be less useful if you don't respect the orders of infinity system 
#### FORMALISATION: Co-germs, co-germ differentials
- Explain what a "germ" is. Curves could agree everywhere it matters but technically be different curves. A "germ" is a collection of curves which should be treated the same way.
- A variable is a co-germ. differentials are co-germ differential 1-forms (REMARK: this is how Toby Bartels defines it. CX => R)
- df works as advertised, two curves with the same germ will give the same result. REMARK: Bartels assumes f is a global function f:X-R. QUERY: do we need to be using different terminology?
- REMARK: Bartels treats ds as a cogerm differential form because you can define any phi(forms of the form df) as a cogerm differential form. I.e. ds is a function of dx and dy so it a cogerm differential form, by lifting in the obvious way. I, on the other hand, get to ds more directly. Though the result is the same (is it always the same?)
- Not 1-forms because they have access to all rates of change of the curve, 1-forms only work with equivalence classes with the same tangent vector.

### Episode: partial differentials
- PLANNING
- Can we do all of the above for partials?
- The ideal would be to have an object "partial y" and another "partial x" so that one can get the partial differential by dividing one by the other. We would have to have a way of "keeping fixed" a variable.
- You can use a wedge product to perform operations on total differentials to find algebraic identities for partial derivatives: https://johncarlosbaez.wordpress.com/2021/09/13/the-cyclic-identity-for-partial-derivatives/
- Gordon Plotkin has an equational theory: https://math.ucr.edu/home/baez/mathematical/ACTUCR/Plotkin_Partial_Differentiation.pdf
- Baez commentary: https://johncarlosbaez.wordpress.com/2020/05/18/a-complete-axiomatisation-of-partial-differentiation/
#### Complaints and ambiguities
- Arnol'd's complaint
- Functional view
- Lots of "reform" ideas (detail put off to later)
#### In the spirit of Leibniz
- Let's try to do someting in the Leibniz spirit as much as we can. We would like partial f to be a standalone thing which we can take a ratio of.
- Suppose we have f is a function of x and y, an expression in f and its differentials probably defines a surface not a curve. This is not the only possibility (can we perhaps always have surfaces if f is just a variable).
- On a surface there are lots of different directions we could differentiate in, unlike with 1D, where all we didn't know what the parameter's rate. Now we don't know which direction that parameter might be taking us.
- Pro tem, let us try defining partial(f, w). w defines lines of constant w. We would like to be able to differentiate across w (why?) but we can't step between them.
- So assume we have a vector field T of curves parametrised by t that runs across them. Now they are linked up.
- partial(f, w)= pick an arbitrary foliation of curves with parameter t that is not colinear to w lines call it gamma_w(t) and pass that to df (which acts on curves). We get df/dt along gamma.
- This obviously depends on the T foliation. Oh no! But it cancels out just as for total differentials.
- partial(f,w)/partial(f,x) has no such dependency and is a ratio of the two.
- Notationally we can lift partial f so that it is a function on arbitrary gammas, then partial f/partial x (restricted to w) = partial(f,w)/partial(f,x)

## LOAD-BEARING
- [ ] Epsilon-style definition of derivative, Lim epsilon->0 of 1/epsilon( f(x+ε) − f(x)) — it or something like it is needed no later than the definition of df, which will involve some kind of limit. We also need it if we want to contrast differentials v derivatives  — status: idea
- [ ] The Limit approach to delta y/delta x (i.e. with delta x tending to zero). This is the school definition which (we hope) is familiar. It is very close to the epsilon-style definition, but does not absolutely require a functional dependence. It is perhaps more direct than defining df. You can see it on a picture. It is useful (because it explains what is going at school level) — status: idea
- [ ] Limits - do we need an aside on limits, what they are. The best definition (I think) is Hintikka style. If I claim a limit, you can challenge by picking  an accuracy (delta) and I can pick (for a sequence) an N after which I am inside that accuracy.
- [ ] Notation: is dy/dx a function of x? How do you indicate applying it —
      dy/dx(1)? dy/dx|_{x=1}? — and why plain dy/dx is fine most of the time
      *because of* the approach being taken (this resolves the friend's
      "what type is y" objection) 

## DIGRESSION
#### Inexact variables
- Example: s, not defined everyewhere on the plane, but each point on a curve can have s defined.
- If we don't want to assume a base point, we have a "torsor".
- d works fine, because we only care about subtraction.
- Example: W (work in thermodynamics) is worse than s because of loops and that it can subtract as well as add. But this is not a problem in a small enough neighbourhood. Co-germ differentials only care about small neighbourhoods.
- REMARK: mentioning torsors sometimes would be useful because primitives of functions aka "indefinite integrals" are also torsors.
### Curvature
- Curvature via second order — hook: "what does the *next* term buy you. This is a way of showing how one can manipulate differentials at higher order to get a nice formula without thinking about functions at all.
- What curvature is in Newtonian tems (intersections of normals)
- Why that is the same thing as rate of change of angle wrt distance along curve
- Geometrically?". It can illustrate how pleasant it is to play with d. You can just apply d to tan theta=dy/dx and the rest follows very easily algebraically 
  - MATH: tan theta = dy / dx
  - d(tan theta) = sec^2 theta d theta; d(dy/dx) = frac(d x d^2 y - d y d ^2 x, d x^2)
  - sec theta = dx / ds
  - d theta = frac(d x d^2 y - d y d ^2 x, d s^2)
  - ds^2 = dx^2 + dy^2
  - d theta = frac(d x d^2 y - d y d ^2 x, dx^2 + dy^2)
  - d theta / d s =
- Useful to point out that theta and s are have order zero.
### Infinitessimals
- hook: how do we have infinitely short sides without having infinitesimal numbers? "safe," makes minimal ontological demands
- Idea of ininitesimaly small lines making curves polygons (old) v infinitesimally small numbers (newer with Leibniz).
- Can introduce free infinitesimals. Does not tell you anything beyond polynomials in them(?). If we had a model we could do more.
-  Hardy on growth-rate infinitesimals (early 20th c. calculus texts defining infinitesimals as classes of function by growth rate) 
- Robinson / nonstandard analysis — hook: still needed/ Why do we need this?
- Nilpotents, and Grothendieck's desire for a "double point" — hook: a genuinely motivating picture for why algebraic geometry wants nilpotents at all
### MISC
- [ ] Why 1-forms don't work (no access to higher derivatives), why forms really don't work (even if integrable, you can make a foliation and then convert all forms to scalars via Hodge duality and then take a d of them, but you need parametrisation for all that - is this true? Can it be fleshed out. Note: it was my original model but it may illustrate why you can't get away from the parameter).

- [ ] A bonus point is that things like ds =dx^2 +dy^2 are natural and in SR you can just do tau = (similar) and derive d tau/dt = gamma (or 1/gamma) directly again just using algebra. All type checks nicely (see inexactnes. tau is a famous example of an inexact quantity. See twin paradox (a future physics blog).

-
- [ ] Once we explore "things like ds" we realise we can do inexact differentials, but we didn't know what they were.
- [ ] Explain exactness simply via thing like potentials and voltage. This leads to small-delta notation
- [ ] Small-delta notation for inexact differentials → brief discussion of exactness — hook: still needed — status: idea
- [ ] Bigger digression: thermodynamics. dU=delta q + delta W (how you can guess an exact differential exists if you now something about an inexact one - in this case W along adiabatic paths)
- [ ] Exactness → homology (resist going further than a taste) — hook: still 0needed, and a place to explicitly signal "you can stop here" — status: idea

## OPEN QUESTIONS (mine, not the reader's)
- Does 1D integration belong as part of the main discussion and then 2+D integration belong later?
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
### Dimensional analysis
Comment: Order of Infinity on the spine will be a sort of dimensional analysis. This chunk considers what dimensional analysis "really" is.    
- How does this work?
- What physicists do: have a series of "base dimensions" (dimensional axes?), eg MLT. Every expression has a "dimension" such as MLT^-1 for momentum. ("In physics, these are your chosen base dimensions (Mass, Length, Time). In algebraic geometry and representation theory, these axes are formally called the characters of the group, or simply the coordinate axes of the parameter space"). For the most part we just care about how many there are.
- Formalising this parametrically, we have a *Weight Lattice* here Z^3 with points that might be (1, 1, -1) for momentum. If continuous we say _weight space_. 
- Final values live in the "Target Field" or "Base Field" or Value Ring): This is the immediate mathematical set where the single output number lands. For standard physics, this is simply the field of Real Numbers (\(\mathbb{R}\)) or Complex Numbers (\(\mathbb{C}\)). If you are dealing with infinitesimals (like in your Hardy field or Leibniz model), this would be a specialized field of fractions or germs.
- An "algebraic torus" of values - sometimes called the "scale group" - the scale factors. Here it might be R^3_+ It is positive because you don't scale by zero (no inverse, so not a group) and negative is no good because units can't be negative. For scaling the units. A change of units is the action of an element of the scale group. Also: Structure Group and scaling torus. Eg (1,2,1) means double the second dimension (in MLT that would be length).
- Definition: "algebraic torus" means the multiplicative group of non-zero numbers (C or R)
-  "variables" like x or v take an element of the scale group and give a value in the target field.
-  Expressions involving variable are also functions of the same kind with operations pulled back from the target field. Eg x * v, would be s in scale group -> x(s)*v(s) 
- The Configuration Space (or Ambient State Space): If you combine all your variables together (e.g., you have 5 variables: \(x_1, x_2, x_3, x_4, x_5\)), the values collectively live in a high-dimensional vector space, such as \(\mathbb{R}^{5}\). This is the "Lab Space" of things we might actually measure in the lab. You don't measure expressions you measure variables.
-  The toric ideal is the set of all polynomials that scale by a uniform factor under any element of the scaling group
-  If I=toric ideal, The Variety Z(I) is the Toric Variety.
-  A toric variety is a "a geometric space that contains this scaling group as an open, dense piece" (?). "A Toric Variety is an algebraic variety (a shape defined by polynomial equations) that contains a coordinate scaling group (a torus) as a dense open set, such that the scaling group can act on the entire shape without breaking its geometric structure."
-- Zariski Open: In this topology, a set is "closed" if it is the boundary carved out by a polynomial equation (like a line or a surface). Therefore, an "open set" is just the entire space minus a few lower-dimensional polynomial lines or cuts.
   A subset \(A\) is dense in \(X\) if every non-empty open set in \(X\) contains at least one point of \(A\).
   A subset \(A\) is dense in \(X\) if the only closed set containing \(A\) is the entire space \(X\) itself.
--  Let's look at your physics example: \(vt - d = 0\).If none of your variables are zero (\(v \neq 0, t \neq 0, d \neq 0\)), you can pick any value for \(v\) and \(t\), and \(d\) is instantly locked in by multiplication (\(d = vt\)).
-- This "nonzero chunk" of the surface behaves exactly like a 2D scaling group—you can multiply \(v\) and \(t\) by any non-zero scale factor and cleanly glide across the surface. This is your algebraic torus.
-- What happens if a variable is zero? If \(t = 0\), then \(d = 0\), and \(v\) can be absolutely anything. This forms a boundary line on the edge of your surface.
-  The entire surface (the Toric Variety) is made of that smooth, non-zero scaling region (the dense open torus) plus the sharp boundary lines where things hit zero.
-  Digression: how to make a toric variety from fans. A fan is a collection of cones drawn on a lattice.

### torsors
Needed or useful for dimensional analysis. Unclear whether we need it.
-  Group that has forgotten its identity element
- Eg, Voltage, Energy, date, notes. You can't add these together, you can subtract them and get a torsor. You can add a torsor to one of these to get another one of these values. We could have used torsors for our values.
- John Baez has a very good informal discussion: https://math.ucr.edu/home/baez/torsors.html

### Doing parametric dimensions for our story
- Base dimension - The pacing of the curve (t). This is a torsor because t1+t2 is meaningless (can't add points)
- Scale group Diff(R) - group of smooth changes of variable t=g(s)
- Target field - R
- Weight Lattice - Z (order of differential)
- Toric Ideal => contact ideal. Formulas like dy -d'(x)dx
- Toric Variety => Jet Variety
- Resulting orbit space (quotienting) is intrinsic geometry.
### Maxwell's thermodynamics
- Could be an interesting illustration of several things.
- Thorough work up of geometry: https://arxiv.org/pdf/1811.04227
- More informal on graphical methods: https://academicweb.nd.edu/~powers/ame.20231/gibbs1873a.pdf?utm_campaign=redirect&utm_medium=web&utm_source=www3.nd.edu
- 

## Reading List
- https://mathoverflow.net/questions/65829/de-rham-cohomology-vs-iterated-tangent-bundles
- P.-A. Meyer, Qu'est ce qu'une différentielle d'ordre  n  , Exposition. Math. 7 (1989), 249–264.
- Laksov, Dan; Thorup, Anders, These are the differentials of order  n  . Trans. Amer. Math. Soc. 351 (1999), no. 4, 1293–1353. https://www.ams.org/journals/tran/1999-351-04/S0002-9947-99-02120-0/S0002-9947-99-02120-0.pdf
- https://ncatlab.org/nlab/show/K%C3%A4hler+differential
- https://math.stackexchange.com/questions/2843918/the-second-differential-versus-the-differential-of-a-differential-form
- "Total and Partial Differentials as Algebraically Manipulable Entities" https://arxiv.org/pdf/2210.07958 This defines d as an operator in much the same style as in this piece. It uses hyperreals as a model. It assumes all variables depend on some "underlying variable" $q$. This $q$ is our parameter for a curve. There is no connection made between this and the study of curves. Otherwise the discussion is much the same. No discussion of integration, but some discussion of partials.
### partial derivatives
- https://math.stackexchange.com/questions/320228/the-notation-for-partial-derivatives
- https://math.stackexchange.com/questions/2200982/notation-for-partial-derivatives
- https://math.stackexchange.com/questions/1091097/partial-derivative-notation
- https://arxiv.org/html/math/9906079v1
- 

## Old Notes
- [ ] Orders of infinity are probably "graded rings" or "graded fields". Or perhaps the homogenous elements of such a graded field. There may be a relationship to Toric Algebras or Tropical Geometry about which I know nothing. A Hady Field captures growth in a similar way but allows the sums of terms of different growth or order. There is a natural valuation with ultrametric inequality (a Krull valuation?) and some interesting structure including a graded algebra core.
- [ ] Dimensional analysis can be done in more than one way. See Terrence Tao.   Physicists handle dimension differently. Each "dimension" is its own vector space (mass would be 1D). Multiplication is tensor product. Changing units, involves rescaling. You have a scale group and this is what Toric XXX studies. But for us we aren't scaling by numbers but by functions (t), that group is the Jet Group, this gives you a Jet Variety not a Toric variety or a Diffiety. Axes are x, y, dx, dy, d2x etc. Cones are the COntact Ideal, dy=f'(x)dx. In a sense physicists dimension are zeroth order and these are dynamic dimensions.

-

