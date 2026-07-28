(() => {
  const script = document.currentScript || Array.from(document.scripts).find((s) => /visitor-counter\.js/.test(s.src));
  const endpoint = new URL("../counter.php", script ? script.src : window.location.href);
  const pagePath = window.location.pathname.replace(/\/index\.html$/i, "/") || "/";

  function formatNumber(value) {
    return new Intl.NumberFormat("en-US").format(Number(value || 0));
  }

  function mountCounter(data) {
    const target = document.querySelector(".site-footer__bottom") || document.querySelector("footer") || document.body;
    if (!target) return;

    let counter = document.querySelector("[data-umtts-counter]");
    if (!counter) {
      counter = document.createElement("div");
      counter.className = "umtts-counter";
      counter.setAttribute("data-umtts-counter", "true");
      target.appendChild(counter);
    }

    counter.innerHTML = `
      <span>Site visitors: ${formatNumber(data.site_visitors)}</span>
      <span>Page visitors: ${formatNumber(data.page_visitors)}</span>
    `;
  }

  fetch(`${endpoint.href}?page=${encodeURIComponent(pagePath)}`, {
    method: "GET",
    credentials: "include",
    cache: "no-store"
  })
    .then((response) => response.ok ? response.json() : Promise.reject(response.status))
    .then(mountCounter)
    .catch(() => {
      const target = document.querySelector(".site-footer__bottom") || document.querySelector("footer");
      if (!target) return;
      const counter = document.createElement("div");
      counter.className = "umtts-counter umtts-counter--offline";
      counter.textContent = "Visitor counter unavailable on this host.";
      target.appendChild(counter);
    });
})();
