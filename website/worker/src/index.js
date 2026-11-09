// Receives the Download page form and adds the visitor to the Mailchimp
// audience. Visitors who ask for release emails are subscribed; everyone else
// is stored as a non-subscribed ("transactional") contact. Existing subscribers
// stay subscribed if they leave the box unticked.
import { createHash } from "node:crypto";

const TAG = "pioneer-download";

export default {
  async fetch(request, env) {
    const origin = request.headers.get("Origin") ?? "";
    const allowed = env.ALLOWED_ORIGINS.split(",").map((o) => o.trim());
    const cors = allowed.includes(origin)
      ? { "Access-Control-Allow-Origin": origin, Vary: "Origin" }
      : {};

    if (request.method === "OPTIONS") {
      return new Response(null, {
        status: 204,
        headers: { ...cors, "Access-Control-Allow-Methods": "POST", "Access-Control-Allow-Headers": "Content-Type" },
      });
    }
    if (request.method !== "POST") return new Response("Method not allowed", { status: 405 });
    if (!cors["Access-Control-Allow-Origin"]) return new Response("Forbidden", { status: 403 });

    let body;
    try {
      body = await request.json();
    } catch {
      return json({ error: "Invalid JSON" }, 400, cors);
    }

    // Bots fill the hidden field; accept silently and do nothing.
    if (body.website) return json({ ok: true }, 200, cors);

    const name = clean(body.name, 200);
    const email = clean(body.email, 254).toLowerCase();
    const organization = clean(body.organization, 200);
    const updates = body.updates === true;
    if (!name || !organization || !/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(email)) {
      return json({ error: "Missing or invalid fields" }, 400, cors);
    }

    const [first, ...rest] = name.split(/\s+/);
    const member = {
      email_address: email,
      status_if_new: updates ? "subscribed" : "transactional",
      merge_fields: { FNAME: first, LNAME: rest.join(" "), ORG: organization },
    };
    if (updates) member.status = "subscribed";

    const hash = createHash("md5").update(email).digest("hex");
    const res = await mailchimp(env, "PUT", `members/${hash}`, member);
    if (!res.ok) {
      const detail = await res.text();
      console.error("Mailchimp member update failed", res.status, detail, { name, email, organization, updates });
      // Mailchimp refuses addresses it considers fake or mistyped; the form
      // asks the visitor to correct it.
      if (res.status === 400 && /email address/i.test(detail)) {
        return json({ error: "invalid_email" }, 422, cors);
      }
      return json({ error: "Could not save" }, 502, cors);
    }
    const tagRes = await mailchimp(env, "POST", `members/${hash}/tags`, {
      tags: [{ name: TAG, status: "active" }],
    });
    if (!tagRes.ok) console.error("Mailchimp tag failed", tagRes.status, await tagRes.text());

    return json({ ok: true }, 200, cors);
  },
};

function clean(value, max) {
  return typeof value === "string" ? value.trim().slice(0, max) : "";
}

function json(data, status, headers) {
  return new Response(JSON.stringify(data), {
    status,
    headers: { ...headers, "Content-Type": "application/json" },
  });
}

function mailchimp(env, method, path, body) {
  const key = env.MAILCHIMP_API_KEY.trim();
  return fetch(`https://${env.MAILCHIMP_DC}.api.mailchimp.com/3.0/lists/${env.MAILCHIMP_LIST_ID}/${path}`, {
    method,
    headers: {
      Authorization: `Basic ${btoa(`pioneer:${key}`)}`,
      "Content-Type": "application/json",
    },
    body: JSON.stringify(body),
  });
}
