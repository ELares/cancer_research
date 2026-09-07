# Preprint full text: what is reachable, and what is not

Cancer preprints are the largest block of literature this project can
identify but mostly cannot read. Counts are Europe PMC's own
`hitCount`, measured by `scripts/preprint_access_check.py`.

| set | records |
|---|--:|
| cancer preprints Europe PMC indexes | 108,400 |
| of those, whose full text Europe PMC holds | 16,436 |
| share | 15.2% |

So roughly 91,964 cancer preprints are identified and unread.

## Where the rest live, and why they stay there

bioRxiv and medRxiv publish a metadata API that names a full-text
location for every record -- a `jatsxml` URL on their own domain. That
endpoint needs no account and is NOT paywalled.

It is also, from this project's client, unavailable: `HTTP 429`
(Cloudflare error 1015). Retried after 20s and again after 45s: the
same. A 429 that does not clear across a minute of backoff is a
standing block on automated access, not a transient rate limit.

**That block is respected.** It would be straightforward to send a
browser User-Agent and get a different answer, and doing so would
circumvent an access control the operator deliberately put in place.
The corpus does not do that anywhere, and this page exists so nobody
has to rediscover the temptation.

## The route that does exist costs money

bioRxiv's supported bulk text-and-data-mining channel is a
requester-pays Amazon S3 bucket: the data is free, the transfer is
billed to whoever asks for it. That is a spending decision, not a
technical one, so it is recorded here rather than taken.

## A correction

An earlier note in this project said bioRxiv full text was available
ONLY through requester-pays S3. That was wrong in a way worth stating:
the per-article JATS XML is genuinely free and public, and the
requester-pays bucket is the BULK route. What blocks the per-article
route is bot protection, not licensing and not cost. The distinction
matters because the two have different remedies -- one needs a
conversation with bioRxiv, the other needs a budget.

## What this means for the corpus

The 16,436 preprints Europe PMC does hold are already
collected by `scripts/corpus_expand_fetch.py`, which asks Europe PMC and
accepts its refusal. The remaining 91,964 are held as metadata and
abstracts, which are themselves worth having: bibliographic data and
abstracts are generally not copyrightable (Feist), and an abstract is
enough for the census questions this project asks most often.
